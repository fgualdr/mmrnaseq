# mm_realize_bam.awk
#
# First input: decision table (resolution.tsv plus EXCLUDED fragment rows).
# Second input: query-name-sorted SAM with header.
#
# Realizes three deliberately distinct BAM semantics:
#   full_out         unique + selected FDR-resolved + every unresolved/excluded alignment
#   quantitative_out unique + every globally expected-FDR-resolved fragment
#   highconf_out     unique + resolved fragments with individual pmax >= highconf_min
#
# IMPORTANT: ZP is emitted only on a selected resolved alignment, where it is the
# posterior probability of that selected candidate. Unresolved alternatives do not
# receive ZP because resolution.tsv contains only fragment pmax, not each alternative's
# candidate posterior; attaching pmax to every alternative would be misleading.

BEGIN {
    FS = OFS = "\t"
    last_frag = ""
    unique_frag = resolved_fdr_frag = resolved_highconf_frag = 0
    resolved_lowerconf_frag = unresolved_frag = excluded_frag = 0
    full_records = quantitative_records = highconf_records = 0
    if (highconf_min == "") highconf_min = 0.95
}

function band(x,y,    r,p,xb,yb) {
    r=0; p=1
    while (x>0 || y>0) {
        xb=x%2; yb=y%2
        if (xb>=1 && yb>=1) r+=p
        x=int(x/2); y=int(y/2); p*=2
    }
    return r
}
function clearbit(x,b) { return x-band(x,b) }
function tagv(prefix,    i,a) {
    for (i=12; i<=NF; i++) {
        if (index($i,prefix)==1) {
            split($i,a,":")
            return a[3]
        }
    }
    return ""
}
function clean_optional_tags(    i,s) {
    s=""
    for (i=12; i<=NF; i++) {
        if ($i ~ /^(NH:i:|HI:i:|ZR:Z:|ZP:f:|ZO:i:)/) continue
        s = s OFS $i
    }
    return s
}
function append_passthrough(out,status_s,orig_nh,    i,s) {
    s=""
    for (i=1; i<=11; i++) s=s (i==1?"":OFS) $i
    for (i=12; i<=NF; i++) {
        if ($i ~ /^(ZR:Z:|ZP:f:|ZO:i:)/) continue
        s=s OFS $i
    }
    s=s OFS "ZR:Z:" status_s OFS "ZO:i:" orig_nh
    print s > out
}
function posterior_mapq(p,    q) {
    if (p >= 1.0) return 60
    if (p <= 0.0) return 0
    q = int(-10.0 * log(1.0-p) / log(10.0) + 0.5)
    if (q > 60) q=60
    if (q < 0) q=0
    return q
}
function append_resolved(out,status_s,p,orig_nh,    i,s,flag,q) {
    # Selected STAR alternatives may carry the secondary bit. Once the Bayesian
    # resolver selects one HI, clear only 0x100 and preserve all pairing/orientation
    # metadata (mate flags, CIGAR, RNEXT/PNEXT/TLEN, duplicate flag, etc.).
    flag=clearbit($2,256)
    q=posterior_mapq(p)
    s=$1 OFS flag OFS $3 OFS $4 OFS q OFS $6 OFS $7 OFS $8 OFS $9 OFS $10 OFS $11
    for (i=12; i<=NF; i++) {
        if ($i ~ /^(NH:i:|HI:i:|ZR:Z:|ZP:f:|ZO:i:)/) continue
        s=s OFS $i
    }
    s=s OFS "NH:i:1" OFS "HI:i:1" OFS "ZR:Z:" status_s \
        OFS "ZP:f:" sprintf("%.8g",p) OFS "ZO:i:" orig_nh
    print s > out
}

NR==FNR {
    if (FNR==1) {
        for (i=1; i<=NF; i++) h[$i]=i
        next
    }
    frag=$h["fragment_id"]
    status[frag]=$h["status"]
    selhi[frag]=$h["selected_hi"]+0
    pmax[frag]=$h["pmax"]+0
    next
}

/^@/ {
    print > full_out
    print > quantitative_out
    print > highconf_out
    next
}

{
    frag=$1
    nh=tagv("NH:i:"); if (nh=="") nh=1; nh+=0
    hi=tagv("HI:i:"); if (hi=="") hi=1; hi+=0
    st=status[frag]
    p=pmax[frag]+0

    # Fragment-level QC: classify from the decision table first. An excluded
    # chimeric/supplementary fragment may still have NH=1 and must not be counted
    # as an ordinary trusted unique fragment.
    if (frag != last_frag) {
        if (st=="EXCLUDED") {
            excluded_frag++
        } else if (st=="RESOLVED") {
            resolved_fdr_frag++
            if (p >= highconf_min) resolved_highconf_frag++
            else resolved_lowerconf_frag++
        } else if (st=="UNRESOLVED") {
            unresolved_frag++
        } else if (nh==1) {
            unique_frag++
        } else {
            excluded_frag++
        }
        last_frag=frag
    }

    if (st=="EXCLUDED") {
        append_passthrough(full_out,"EXCLUDED",nh)
        full_records++
        next
    }

    if (nh==1 && st=="") {
        # Trusted unique fragments are unchanged and are present in every tier.
        print > full_out
        print > quantitative_out
        print > highconf_out
        full_records++
        quantitative_records++
        highconf_records++
        next
    }

    if (st=="RESOLVED") {
        if (hi==selhi[frag]) {
            append_resolved(full_out,"RESOLVED_FDR",p,nh)
            append_resolved(quantitative_out,"RESOLVED_FDR",p,nh)
            full_records++
            quantitative_records++
            if (p >= highconf_min) {
                append_resolved(highconf_out,"RESOLVED_HIGHCONF",p,nh)
                highconf_records++
            }
        }
        next
    }

    if (st=="UNRESOLVED") {
        # Preserve every STAR alternative only in the full ambiguity BAM. Do not
        # emit ZP: resolution.tsv supplies fragment pmax rather than this specific
        # alternative's posterior probability.
        append_passthrough(full_out,"UNRESOLVED",nh)
        full_records++
        next
    }

    # Any multimapping record without a valid decision is conservatively retained
    # only in the archival full BAM and labelled EXCLUDED.
    append_passthrough(full_out,"EXCLUDED",nh)
    full_records++
}

END {
    print "metric","value" > qc_out
    print "unique_fragments",unique_frag > qc_out
    print "resolved_fdr_fragments",resolved_fdr_frag > qc_out
    print "resolved_highconf_fragments",resolved_highconf_frag > qc_out
    print "resolved_lower_confidence_fragments",resolved_lowerconf_frag > qc_out
    print "unresolved_fragments",unresolved_frag > qc_out
    print "excluded_fragments",excluded_frag > qc_out
    print "highconf_min_pmax",sprintf("%.12g",highconf_min) > qc_out
    print "full_bam_alignment_records",full_records > qc_out
    print "quantitative_bam_alignment_records",quantitative_records > qc_out
    print "highconf_bam_alignment_records",highconf_records > qc_out
}
