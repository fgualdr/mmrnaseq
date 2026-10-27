# mm_build_candidates.awk
# Build exact STAR fragment-candidate and unique-anchor contracts from a name-sorted SAM stream.
# Portable awk: no GNU-only bitwise functions are required.

BEGIN {
    FS = OFS = "\t"
    if (single_end == "") single_end = 0
    if (strandedness == "") strandedness = "unstranded"
    current_q = ""
    total_fragments = unique_fragments = multi_fragments = excluded_fragments = 0
    candidate_rows = anchor_rows = 0
    inconsistent_pair_tags = invalid_pair_orientation = incomplete_pair = 0

    print "fragment_id","candidate_id","hi","nh","paired","orientation_group","as","nM","nm","query_bases","aligned_query_bases","normalized_as","r1_chrom","r1_pos1","r1_flag","r1_cigar","r1_blocks","r2_chrom","r2_pos1","r2_flag","r2_cigar","r2_blocks","candidate_signature" > candidate_out
    print "fragment_id","orientation_group","as","nM","nm","query_bases","aligned_query_bases","normalized_as","paired" > unique_out
    print "fragment_id","chrom","pos0","orientation_group","weight","mate","block_index","query_bases" > anchor_out
    print "fragment_id","reason" > excluded_out
}

function band(x, y,    r,p,xb,yb) {
    r=0; p=1
    while (x>0 || y>0) {
        xb=x%2; yb=y%2
        if (xb>=1 && yb>=1) r+=p
        x=int(x/2); y=int(y/2); p*=2
    }
    return r
}
function hasflag(flag, bit) { return band(flag, bit) != 0 }
function absval(x) { return x < 0 ? -x : x }
function tag_value(prefix,    i,a) {
    for (i=12; i<=NF; i++) {
        if (index($i, prefix)==1) {
            split($i,a,":")
            return a[3]
        }
    }
    return ""
}
function has_tag(prefix,    i) {
    for (i=12; i<=NF; i++) if (index($i,prefix)==1) return 1
    return 0
}
function cigar_query_length(cigar,    s,tok,n,op,total) {
    if (cigar=="*" || cigar=="") return 0
    s=cigar; total=0
    while (match(s,/[0-9]+[MIDNSHP=X]/)) {
        tok=substr(s,RSTART,RLENGTH); n=substr(tok,1,length(tok)-1)+0; op=substr(tok,length(tok),1)
        if (op=="M" || op=="I" || op=="S" || op=="=" || op=="X") total+=n
        s=substr(s,RSTART+RLENGTH)
    }
    return total
}
function cigar_aligned_query(cigar,    s,tok,n,op,total) {
    if (cigar=="*" || cigar=="") return 0
    s=cigar; total=0
    while (match(s,/[0-9]+[MIDNSHP=X]/)) {
        tok=substr(s,RSTART,RLENGTH); n=substr(tok,1,length(tok)-1)+0; op=substr(tok,length(tok),1)
        if (op=="M" || op=="I" || op=="=" || op=="X") total+=n
        s=substr(s,RSTART+RLENGTH)
    }
    return total
}
function cigar_nm_gaps(cigar,    s,tok,n,op,go,gl) {
    s=cigar; go=0; gl=0
    while (match(s,/[0-9]+[MIDNSHP=X]/)) {
        tok=substr(s,RSTART,RLENGTH); n=substr(tok,1,length(tok)-1)+0; op=substr(tok,length(tok),1)
        if (op=="I" || op=="D") {go++; gl+=n}
        s=substr(s,RSTART+RLENGTH)
    }
    return go SUBSEP gl
}
function blocks_string(chrom,pos1,cigar,mass,    total,s,tok,n,op,r,bs,bq,out,idx,w) {
    total=cigar_aligned_query(cigar)
    if (total<=0 || chrom=="*" || pos1<=0) return "."
    s=cigar; r=pos1-1; bs=r; bq=0; out=""; idx=0
    while (match(s,/[0-9]+[MIDNSHP=X]/)) {
        tok=substr(s,RSTART,RLENGTH); n=substr(tok,1,length(tok)-1)+0; op=substr(tok,length(tok),1)
        if (op=="M" || op=="=" || op=="X") {r+=n; bq+=n}
        else if (op=="I") {bq+=n}
        else if (op=="D") {r+=n}
        else if (op=="N") {
            if (bq>0) {
                idx++; w=mass*bq/total
                if (out!="") out=out"," 
                out=out chrom ":" bs "-" r ":" sprintf("%.12g",w)
            }
            r+=n; bs=r; bq=0
        }
        s=substr(s,RSTART+RLENGTH)
    }
    if (bq>0) {
        idx++; w=mass*bq/total
        if (out!="") out=out"," 
        out=out chrom ":" bs "-" r ":" sprintf("%.12g",w)
    }
    return out=="" ? "." : out
}
function emit_anchor_blocks(fragment,blocks,orientation,mate,query_bases,    n,a,i,b,c,se,sw,mid) {
    if (blocks=="." || blocks=="") return
    n=split(blocks,a,",")
    for (i=1;i<=n;i++) {
        split(a[i],b,":"); split(b[2],c,"-")
        mid=int(((c[1]+0)+(c[2]+0)-1)/2)
        print fragment,b[1],mid,orientation,b[3],mate,i,query_bases > anchor_out
        anchor_rows++
    }
}
function observed_orientation(flag1,flag2,paired,    g1,g2) {
    if (strandedness=="unstranded" || strandedness=="non_specific") return "."
    if (!paired) return hasflag(flag1,16) ? "rw" : "fw"
    g1=hasflag(flag1,16) ? "fw" : "rw"
    g2=hasflag(flag2,16) ? "rw" : "fw"
    if (g1!=g2) return "INVALID"
    return g1
}
function clear_group(    k) {
    delete rflag; delete rchrom; delete rpos; delete rcigar; delete rseq
    delete ras; delete rnm; delete rNM; delete rnh; delete present
    delete his; maxhi=0; fragment_supp=0; fragment_sa=0
}
function store_record(    hi,nh,as,nm,NM,mate,key) {
    hi=tag_value("HI:i:"); if (hi=="") hi=1; hi+=0
    nh=tag_value("NH:i:"); if (nh=="") nh=1; nh+=0
    as=tag_value("AS:i:"); if (as=="") as=0; as+=0
    nm=tag_value("nM:i:"); if (nm=="") nm=tag_value("NM:i:"); if (nm=="") nm=0; nm+=0
    NM=tag_value("NM:i:"); if (NM=="") NM=nm; NM+=0
    mate=0
    if (!single_end) {
        if (hasflag($2,64)) mate=1
        else if (hasflag($2,128)) mate=2
        else mate=9
    }
    if (hasflag($2,2048)) fragment_supp=1
    if (has_tag("SA:Z:")) fragment_sa=1
    key=hi SUBSEP mate
    if (!present[key] || !hasflag($2,2048)) {
        present[key]=1; rflag[key]=$2+0; rchrom[key]=$3; rpos[key]=$4+0; rcigar[key]=$6; rseq[key]=$10
        ras[key]=as; rnm[key]=nm; rNM[key]=NM; rnh[key]=nh
    }
    his[hi]=1; if (hi>maxhi) maxhi=hi
}
function flush_fragment(    hi,k1,k2,paired,nh,as,nm,NM,q1,q2,aq1,aq2,qtot,aqtot,ori,b1,b2,norm,sig,cid,valid,reason,mass1,mass2) {
    if (current_q=="") return
    total_fragments++
    paired = single_end ? 0 : 1
    if (fragment_supp || fragment_sa) {
        reason = fragment_supp ? "supplementary_alignment" : "chimeric_SA_tag"
        print current_q,reason > excluded_out; excluded_fragments++; clear_group(); return
    }
    valid=0
    for (hi=1;hi<=maxhi;hi++) if (his[hi]) {
        if (single_end) {
            k1=hi SUBSEP 0
            if (!present[k1]) continue
            nh=rnh[k1]; as=ras[k1]; nm=rnm[k1]; NM=rNM[k1]
            q1=(rseq[k1]!="*" ? length(rseq[k1]) : cigar_query_length(rcigar[k1])); q2=0
            aq1=cigar_aligned_query(rcigar[k1]); aq2=0; qtot=q1; aqtot=aq1
            ori=observed_orientation(rflag[k1],0,0)
            b1=blocks_string(rchrom[k1],rpos[k1],rcigar[k1],1.0); b2="."
            mass1=1.0; mass2=0
        } else {
            k1=hi SUBSEP 1; k2=hi SUBSEP 2
            if (!present[k1] || !present[k2]) {incomplete_pair++; continue}
            if (rnh[k1]!=rnh[k2] || ras[k1]!=ras[k2] || rnm[k1]!=rnm[k2]) inconsistent_pair_tags++
            nh=(rnh[k1]>rnh[k2]?rnh[k1]:rnh[k2]); as=(ras[k1]>ras[k2]?ras[k1]:ras[k2]); nm=(rnm[k1]>rnm[k2]?rnm[k1]:rnm[k2]); NM=rNM[k1]+rNM[k2]
            q1=(rseq[k1]!="*" ? length(rseq[k1]) : cigar_query_length(rcigar[k1])); q2=(rseq[k2]!="*" ? length(rseq[k2]) : cigar_query_length(rcigar[k2]))
            aq1=cigar_aligned_query(rcigar[k1]); aq2=cigar_aligned_query(rcigar[k2]); qtot=q1+q2; aqtot=aq1+aq2
            ori=observed_orientation(rflag[k1],rflag[k2],1)
            if (ori=="INVALID") {invalid_pair_orientation++; continue}
            b1=blocks_string(rchrom[k1],rpos[k1],rcigar[k1],0.5); b2=blocks_string(rchrom[k2],rpos[k2],rcigar[k2],0.5)
            mass1=0.5; mass2=0.5
        }
        norm=(qtot>0 ? as/qtot : 0)
        sig=rchrom[k1]":"rpos[k1]":"rflag[k1]":"rcigar[k1]
        if (paired) sig=sig"|"rchrom[k2]":"rpos[k2]":"rflag[k2]":"rcigar[k2]
        cid=current_q"|HI="hi
        if (nh>1) {
            print current_q,cid,hi,nh,paired,ori,as,nm,NM,qtot,aqtot,sprintf("%.12g",norm),rchrom[k1],rpos[k1],rflag[k1],rcigar[k1],b1,(paired?rchrom[k2]:"."),(paired?rpos[k2]:0),(paired?rflag[k2]:0),(paired?rcigar[k2]:"."),b2,sig > candidate_out
            candidate_rows++; valid=1
        } else if (hi==1) {
            print current_q,ori,as,nm,NM,qtot,aqtot,sprintf("%.12g",norm),paired > unique_out
            emit_anchor_blocks(current_q,b1,ori,(paired?"R1":"SE"),qtot)
            if (paired) emit_anchor_blocks(current_q,b2,ori,"R2",qtot)
            valid=1
        }
    }
    # Count fragment once according to NH from any retained record.
    nh=0
    for (hi=1;hi<=maxhi;hi++) if (his[hi]) {
        k1=hi SUBSEP (single_end?0:1)
        if (present[k1]) {nh=rnh[k1]; break}
    }
    if (valid) {
        if (nh>1) multi_fragments++; else unique_fragments++
    } else {
        print current_q,"no_complete_eligible_candidate" > excluded_out; excluded_fragments++
    }
    clear_group()
}

/^@/ { next }
{
    if (current_q!="" && $1!=current_q) flush_fragment()
    if ($1!=current_q) { current_q=$1; clear_group() }
    if (!hasflag($2,4)) store_record()
}
END {
    flush_fragment()
    print "metric","value" > qc_out
    print "total_fragments",total_fragments > qc_out
    print "unique_fragments",unique_fragments > qc_out
    print "multimapper_fragments",multi_fragments > qc_out
    print "excluded_fragments",excluded_fragments > qc_out
    print "candidate_rows",candidate_rows > qc_out
    print "anchor_rows",anchor_rows > qc_out
    print "inconsistent_pair_tags",inconsistent_pair_tags > qc_out
    print "invalid_pair_orientation",invalid_pair_orientation > qc_out
    print "incomplete_pair_candidates",incomplete_pair > qc_out
}
