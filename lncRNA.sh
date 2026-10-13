#!/usr/bin/env bash
# de novo lncRNA discovery from a StringTie GTF (shell/awk re-write of PLAR)
#
# Requirements: bash, awk, bedtools (>=2.26), CPC2, HMMER3 (hmmsearch) + Pfam-A

set -euo pipefail

VERSION="2.0.0"
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

########################## configuration ##########################
# Paths can be edited here, set via environment variables, or given as options.
BED=${LNC_BED_DIR:-$SCRIPT_DIR/bed}
genome=${LNC_GENOME:-/home/quanyi/genome/hg19/GRCh37.p13.genome.fa}
CPC2=${LNC_CPC2:-/home/quanyi/app/CPC2-beta/bin/CPC2.py}
Pfam=${LNC_PFAM:-/home/quanyi/app/PLAR/Pfam/Pfam-A.hmm}

s_fpkm=1          # single-exon FPKM cutoff
m_fpkm=0.1        # multi-exon FPKM cutoff
r_fpkm=0.1        # known (reference) lncRNA FPKM cutoff
min_len=200       # minimal transcript length (sum of exons)
repeat_frac=0.5   # minimal fraction of an exon covered by one repeat/gap feature
s_dist=2000       # minimal distance of single-exon transcripts to protein-coding genes
m_dist=500        # minimal distance of multi-exon transcripts to protein-coding genes
evalue=0.001      # HMMER per-sequence E-value cutoff
threads=$(nproc 2>/dev/null || echo 4); [ "$threads" -gt 16 ] && threads=16
outdir=.
known_list=
mod=
keep_tmp=0
###################################################################

help(){
    cat <<-EOF
	Usage: lncRNA.sh <options> -g|-n <transcripts.gtf>

	### INPUT: GTF file (stringtie output) ###
	### bedtools/CPC2/HMMER(Pfam-A) required ###

	Options:
	    -g  GENCODE reference was used for transcript assembly
	    -n  NONCODE reference was used for transcript assembly
	    -s  single-exon FPKM cutoff (default: $s_fpkm)
	    -m  multi-exon FPKM cutoff (default: $m_fpkm)
	    -r  referenced transcript FPKM cutoff (default: $r_fpkm)
	    -l  known lncRNA transcript ID list, plain or .gz
	        (default: bed/gencode.v29lift37.long_noncoding_RNAs.list.gz for -g,
	         bed/NONCODE.list.gz for -n)
	    -o  output directory (default: ./)
	    -t  threads for CPC2/HMMER (default: $threads)
	    -G  genome FASTA (default: \$LNC_GENOME or $genome)
	    -C  CPC2.py command (default: \$LNC_CPC2 or CPC2.py in PATH)
	    -P  Pfam-A.hmm (default: \$LNC_PFAM or $Pfam)
	    -k  keep intermediate files
	    -v  print version
	    -h  print this help message
	EOF
}

die(){ echo "ERROR: $*" >&2; exit 1; }

# no ARGs error
if [ $# -lt 1 ];then
    help
    exit 1
fi

while getopts "hvgnks:m:r:l:o:t:G:C:P:" arg
do
    case $arg in
        s) s_fpkm=$OPTARG;;
        m) m_fpkm=$OPTARG;;
        r) r_fpkm=$OPTARG;;
        l) known_list=$OPTARG;;
        o) outdir=$OPTARG;;
        t) threads=$OPTARG;;
        G) genome=$OPTARG;;
        C) CPC2=$OPTARG;;
        P) Pfam=$OPTARG;;
        k) keep_tmp=1;;
        g) mod='gencode';;
        n) mod='noncode';;
        v) echo "lncRNA.sh $VERSION"; exit 0;;
        h) help; exit 0;;
        ?) help; exit 1;;
    esac
done
shift $((OPTIND - 1))

[ $# -eq 1 ] || { help; exit 1; }
input=$1
[ -n "$mod" ] || die "Reference has to be either GENCODE (-g) or NONCODE (-n)!"
[ -s "$input" ] || die "input GTF '$input' not found or empty"
for v in "$s_fpkm" "$m_fpkm" "$r_fpkm"; do
    [[ $v =~ ^[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$ ]] || die "FPKM cutoff '$v' is not a number"
done
[[ $threads =~ ^[1-9][0-9]*$ ]] || die "threads '$threads' is not a positive integer"

prefix=$(basename "$input" .gtf)

# BED resources may be gzipped (as shipped) or not (gunzipped by older versions)
resource(){ # $1=file name without .gz
    if [ -s "$BED/$1" ]; then echo "$BED/$1"
    elif [ -s "$BED/$1.gz" ]; then echo "$BED/$1.gz"
    else die "resource $BED/$1(.gz) not found"; fi
}
read_any(){ case $1 in *.gz) gzip -dc -- "$1";; *) cat -- "$1";; esac; }

if [ "$mod" = gencode ]; then
    keyw='ENST'; : "${known_list:=$(resource gencode.v29lift37.long_noncoding_RNAs.list)}"
else
    keyw='NONHSAT'; : "${known_list:=$(resource NONCODE.list)}"
fi
[ -s "$known_list" ] || die "known lncRNA list '$known_list' not found"
repeat_bed=$(resource hg19.repeat_gap.bed)
cds_bed=$(resource CDS.bed)
gene_bed=$(resource hg19_pc_txSE.bed)

# check dependencies before doing any work
command -v bedtools >/dev/null || die "bedtools not found in PATH"
command -v hmmsearch >/dev/null || die "hmmsearch (HMMER3) not found in PATH"
if [ ! -e "${CPC2%% *}" ] && ! command -v "${CPC2%% *}" >/dev/null; then
    command -v CPC2.py >/dev/null || die "CPC2 not found: '$CPC2' (set -C or \$LNC_CPC2)"
    CPC2=$(command -v CPC2.py)
fi
[ -s "$genome" ] || die "genome FASTA '$genome' not found (set -G or \$LNC_GENOME)"
[ -s "$Pfam" ] || die "Pfam-A HMM '$Pfam' not found (set -P or \$LNC_PFAM)"

mkdir -p "$outdir"
outdir=$(cd "$outdir" && pwd)
input=$(cd "$(dirname "$input")" && pwd)/$(basename "$input")
tmp=$(mktemp -d "$outdir/.${prefix}.lncRNA.XXXXXX")
cleanup(){
    local rc=$?
    if [ $rc -ne 0 ]; then
        echo "Run failed, intermediate files kept in $tmp" >&2
    elif [ $keep_tmp -eq 0 ]; then
        rm -rf "$tmp"
    fi
}
trap cleanup EXIT
export LC_ALL=C

# a and b are BED6 + BED6: hit when either strand is unknown or both are equal
STRAND_OK='$6=="." || $12=="." || $6==$12'

###################################################################
# 1. parse GTF, classify transcripts, filter 1 (FPKM and length)
#    tx.tsv: id chr start end strand type n_exon length FPKM source status
###################################################################
parse_gtf(){
    read_any "$known_list" > "$tmp/known.list"
    awk -F'\t' -v OFS='\t' -v mod="$mod" -v keyw="$keyw" \
        -v sf="$s_fpkm" -v mf="$m_fpkm" -v rf="$r_fpkm" -v minlen="$min_len" \
        -v tx="$tmp/tx.tsv" -v exon="$tmp/exon.bed" '
    function attr(a, k,   p, r){
        p = index(a, " " k " \""); if (!p) return ""
        r = substr(a, p + length(k) + 3)
        return substr(r, 1, index(r, "\"") - 1)
    }
    FILENAME == ARGV[1] { if ($1 !~ /^#/ && !($1 in known)) { known[$1]; nknown++ }; next }
    /^#/ || ($3 != "transcript" && $3 != "exon") { next }
    {
        a = " " $9; id = attr(a, "transcript_id")
        if (id == "") next
        if ($3 == "transcript") {
            if (!(id in chr)) ids[++n] = id
            chr[id] = $1; st[id] = $4 - 1; en[id] = $5; sd[id] = $7
            fpkm[id] = attr(a, "FPKM"); ref[id] = attr(a, "reference_id")
            if (fpkm[id] != "") nf++
        } else {
            ne[id]++; len[id] += $5 - $4 + 1
            print $1, $4 - 1, $5, id, 0, $7 > exon
        }
    }
    END{
        if (nf == 0) { print "ERROR: no FPKM attribute found, is this a StringTie GTF?" > "/dev/stderr"; exit 2 }
        for (i = 1; i <= n; i++) {
            id = ids[i]; r = ref[id]; v = fpkm[id] + 0
            if (r == "" && index(id, keyw) == 1) r = id
            type = ne[id] >= 2 ? "multi" : "single"
            status = ""
            if (r != "" && (r in known)) {
                src = "known"; nk++
                if (v <= rf + 0) status = "known_lowFPKM"
                else { nk1++; if (mod == "gencode") status = "known_pass" }
                # known NONCODE transcripts go through the de novo filters as well
            } else if (r != "") {
                src = "reference"; status = "reference_not_lncRNA"
            } else {
                src = "denovo"; nd++; nt[type]++
            }
            if (status == "") {
                nc[type]++
                if (v <= (type == "single" ? sf + 0 : mf + 0) || len[id] < minlen + 0) status = "f1_fpkm_length"
                else if (chr[id] !~ /^chr/) status = "f1_non_chr"
                else np[type]++
            }
            print id, chr[id], st[id], en[id], sd[id], type, ne[id] + 0, len[id] + 0, fpkm[id], src, status > tx
        }
        printf "Total %d transcripts in the input GTF\n", n
        printf "Load %d known lncRNA transcripts from %s\n", nknown, toupper(mod)
        printf "Found %d known lncRNA transcripts in the input GTF\n", nk
        printf "%d out of %d known lncRNA transcripts have low (FPKM<=%s) expression\n", nk - nk1, nk, rf
        printf "Total %d unannotated transcripts\n", nd
        printf "Of these, %d are single-exon and %d are multi-exon transcripts\n", nt["single"], nt["multi"]
        if (mod == "noncode") printf "%d known NONCODE transcripts go through the same filters\n", nc["single"] + nc["multi"] - nd
        printf "%d out of %d single-exon transcripts have FPKM>%s and length>=%dnt\n", np["single"], nc["single"], sf, minlen
        printf "%d out of %d multi-exon transcripts have FPKM>%s and length>=%dnt\n", np["multi"], nc["multi"], mf, minlen
    }' "$tmp/known.list" "$input"
}

###################################################################
# 2-4. repeat/gap, CDS and nearby gene filters, run in parallel
###################################################################
overlap_filters(){
    # candidate transcripts (BED6, score column = type) and their exons
    awk -F'\t' -v OFS='\t' '$11==""{print $2,$3,$4,$1,$6,$5}' "$tmp/tx.tsv" > "$tmp/cand.bed"
    awk -F'\t' -v OFS='\t' 'FILENAME==ARGV[1]{c[$4]; next} $4 in c' "$tmp/cand.bed" "$tmp/exon.bed" > "$tmp/cand_exon.bed"

    # filter 2: any exon covered >= repeat_frac by one repeat/gap feature
    bedtools intersect -a "$tmp/cand_exon.bed" -b "$repeat_bed" -wa -wb -f "$repeat_frac" \
        | awk -F'\t' "$STRAND_OK"'{print $4}' | sort -u > "$tmp/hit.repeat" &
    local p2=$!
    # filter 3: any exon overlapping CDS
    bedtools intersect -a "$tmp/cand_exon.bed" -b "$cds_bed" -wa -wb \
        | awk -F'\t' "$STRAND_OK"'{print $4}' | sort -u > "$tmp/hit.cds" &
    local p3=$!
    # filter 4: protein-coding gene closer than s_dist/m_dist (bedtools closest -d distance)
    local w=$((s_dist > m_dist ? s_dist : m_dist))
    bedtools window -a "$tmp/cand.bed" -b "$gene_bed" -w "$w" \
        | awk -F'\t' -v sd="$s_dist" -v md="$m_dist" "$STRAND_OK"'{
            gap = ($8 > $2 ? $8 : $2) - ($9 < $3 ? $9 : $3)
            if ((gap < 0 ? 0 : gap + 1) < ($5 == "single" ? sd : md)) print $4
        }' | sort -u > "$tmp/hit.gene" &
    local p4=$!
    wait $p2 || die "repeat/gap filter failed"
    wait $p3 || die "CDS filter failed"
    wait $p4 || die "gene distance filter failed"

    update_status f2_repeat "$tmp/hit.repeat" "have overlaps with repeat sequences"
    update_status f3_cds "$tmp/hit.cds" "have overlaps with CDS"
    update_status f4_near_gene "$tmp/hit.gene" "are too close (${s_dist}/${m_dist}bp) to protein-coding genes"
}

update_status(){ # $1=status $2=hit list $3=message; only candidates still passing get updated
    awk -F'\t' -v OFS='\t' -v s="$1" -v msg="$3" -v out="$tmp/tx.tsv.new" '
        FILENAME==ARGV[1]{h[$1]; next}
        $11==""{ n[$6]++; if ($1 in h) { $11 = s; k[$6]++ } }
        { print > out }
        END{
            printf "%d out of %d single-exon transcripts %s\n", k["single"], n["single"], msg
            printf "%d out of %d multi-exon transcripts %s\n", k["multi"], n["multi"], msg
        }' "$2" "$tmp/tx.tsv"
    mv "$tmp/tx.tsv.new" "$tmp/tx.tsv"
}

###################################################################
# 5. coding potential: CPC2 on transcript sequences, then HMMER
#    (3-frame translation vs Pfam-A) on the CPC2 non-coding ones.
#    Transcripts with unknown strand are tested on both strands.
###################################################################
get_fasta(){
    # BED12 of candidate transcripts; unstranded ones twice (+ and - as <id>__rc)
    awk -F'\t' -v OFS='\t' 'FILENAME==ARGV[1]{if ($11=="") c[$1]; next} $4 in c' \
        "$tmp/tx.tsv" "$tmp/cand_exon.bed" | sort -t "$(printf '\t')" -k4,4 -k2,2n \
        | awk -F'\t' -v OFS='\t' '
        function out(   i, sz, bs, s){
            for (i = 1; i <= n; i++) { sz = sz (e[i] - b[i]) ","; bs = bs (b[i] - b[1]) "," }
            s = sd == "." ? "+" : sd
            print c, b[1], end, id, 0, s, b[1], end, 0, n, sz, bs
            if (sd == ".") print c, b[1], end, id "__rc", 0, "-", b[1], end, 0, n, sz, bs
        }
        $4 != id { if (n) out(); id = $4; c = $1; sd = $6; n = 0; end = 0 }
        { b[++n] = $2; e[n] = $3; if ($3 > end) end = $3 }
        END{ if (n) out() }' > "$tmp/cand.bed12"
    [ -s "$tmp/cand.bed12" ] || { : > "$tmp/cand.fa"; return; }
    bedtools getfasta -fi "$genome" -bed "$tmp/cand.bed12" -s -split -name 2>"$tmp/getfasta.log" \
        | sed -e '/^>/s/::.*$//' -e '/^>/s/([+.-])$//' > "$tmp/cand.fa" \
        || die "bedtools getfasta failed, see $tmp/getfasta.log"
    [ -s "$tmp/cand.fa" ] || die "no sequence extracted from $genome (chromosome names differ?)"
}

run_cpc2(){
    local n chunks
    n=$(grep -c '^>' "$tmp/cand.fa" || true)
    chunks=$(( (n + 99) / 100 )); [ $chunks -gt "$threads" ] && chunks=$threads
    mkdir -p "$tmp/cpc2"
    awk -v n="$chunks" -v d="$tmp/cpc2" '/^>/{f = sprintf("%s/chunk_%03d.fa", d, i++ % n)} {print > f}' "$tmp/cand.fa"
    echo "Running CPC2 on $n sequences ($chunks parallel jobs):"
    export LNC_CPC2_CMD=$CPC2
    printf '%s\n' "$tmp"/cpc2/chunk_*.fa | xargs -n 1 -P "$threads" bash -c '
        f=${1%.fa}; $LNC_CPC2_CMD -i "$1" -o "$f" > "$f.log" 2>&1 || { cat "$f.log" >&2; exit 255; }
        [ -s "$f.txt" ] && mv "$f.txt" "$f.out" || mv "$f" "$f.out"' _ \
        || die "CPC2 failed, see $tmp/cpc2/*.log"
    # label column located by header name (CPC2 beta and 1.x differ in extra columns)
    cat "$tmp"/cpc2/chunk_*.out | awk '
        /^#/{ for (i = 1; i <= NF; i++) if ($i == "label") col = i; next }
        { print $1 "\t" $(col ? col : 8) }' > "$tmp/cpc2.tsv"
}

translate3(){ # 3-frame translation of the + strand, standard genetic code
    awk '
    BEGIN{
        aa = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"; nt = "TCAG"
        for (i = 0; i < 64; i++)
            code[substr(nt, int(i / 16) + 1, 1) substr(nt, int(i / 4) % 4 + 1, 1) substr(nt, i % 4 + 1, 1)] = substr(aa, i + 1, 1)
    }
    function out(   f, i, p, c){
        s = toupper(s)
        for (f = 1; f <= 3; f++) {
            p = ""
            for (i = f; i + 2 <= length(s); i += 3) { c = substr(s, i, 3); p = p ((c in code) ? code[c] : "X") }
            if (p != "") print ">" id "|" f "\n" p
        }
    }
    /^>/{ if (id != "") out(); id = substr($1, 2); s = ""; next }
    { s = s $0 }
    END{ if (id != "") out() }' "$@"
}

coding_filters(){
    get_fasta
    if [ ! -s "$tmp/cand.fa" ]; then
        echo "No transcripts left for coding potential filters"
        return
    fi
    run_cpc2
    # a transcript is non-coding only if all its sequences (both strands if unstranded) are
    awk -F'\t' '{ id = $1; sub(/__rc$/, "", id); seen[id]; if ($2 != "noncoding") cod[id] }
        END{ for (id in seen) print id "\t" ((id in cod) ? "coding" : "noncoding") }' "$tmp/cpc2.tsv" > "$tmp/cpc2.call"
    awk -F'\t' '$2=="coding"{print $1}' "$tmp/cpc2.call" > "$tmp/hit.cpc2"
    awk -F'\t' 'FILENAME==ARGV[1]{s[$1]; next} $11=="" && !($1 in s){print $1}' "$tmp/cpc2.call" "$tmp/tx.tsv" > "$tmp/hit.noseq"
    [ -s "$tmp/hit.noseq" ] && echo "WARNING: no genomic sequence for $(wc -l < "$tmp/hit.noseq") transcripts (chromosome names differ from $genome?)"
    update_status f5_no_sequence "$tmp/hit.noseq" "have no genomic sequence" >/dev/null
    update_status f5_cpc2 "$tmp/hit.cpc2" "are coding by CPC2"

    awk -F'\t' 'FILENAME==ARGV[1]{if ($11=="") c[$1]; next}
        /^>/{ id = substr($1, 2); sub(/__rc$/, "", id); p = (id in c) } p' "$tmp/tx.tsv" "$tmp/cand.fa" \
        | translate3 > "$tmp/cand.aa"
    : > "$tmp/hit.pfam"
    [ -s "$tmp/cand.aa" ] && run_hmmer
    update_status f5_pfam "$tmp/hit.pfam" "are coding by HMMER"
}

run_hmmer(){
    # hmmsearch with Z = number of Pfam models reports the same hits and E-values as
    # hmmscan, but is several times faster. hmmsearch threads over target sequence
    # blocks, which barely helps for a few thousand peptides, so the peptides are
    # split into chunks searched by parallel single-threaded processes instead.
    local n z chunks
    n=$(grep -c '^>' "$tmp/cand.aa")
    z=$(grep -c '^NAME ' "$Pfam") || die "no HMM in $Pfam"
    chunks=$(( (n + 99) / 100 )); [ $chunks -gt "$threads" ] && chunks=$threads
    mkdir -p "$tmp/hmmer"
    awk -v n="$chunks" -v d="$tmp/hmmer" '/^>/{f = sprintf("%s/chunk_%03d.aa", d, i++ % n)} {print > f}' "$tmp/cand.aa"
    echo "Running HMMER on $n peptides against $z Pfam-A models ($chunks parallel jobs):"
    export LNC_PFAM_DB=$Pfam LNC_EVALUE=$evalue LNC_Z=$z
    printf '%s\n' "$tmp"/hmmer/chunk_*.aa | xargs -n 1 -P "$threads" bash -c '
        hmmsearch --cpu 1 --noali -o "${1%.aa}.log" -E "$LNC_EVALUE" -Z "$LNC_Z" --domZ "$LNC_Z" \
            --tblout "${1%.aa}.tbl" "$LNC_PFAM_DB" "$1" || { cat "${1%.aa}.log" >&2; exit 255; }' _ \
        || die "hmmsearch failed, see $tmp/hmmer/*.log"
    cat "$tmp"/hmmer/chunk_*.tbl | awk '!/^#/{ id = $1; sub(/\|[123]$/, "", id); sub(/__rc$/, "", id); print id }' \
        | sort -u > "$tmp/hit.pfam"
}

###################################################################
# output
###################################################################
write_outputs(){
    local nadd="$outdir/non_add" f
    local outs=("$outdir/${prefix}_known_lncRNA_f1.gtf" "$outdir/${prefix}_single_exon_f5.gtf" "$outdir/${prefix}_multi_exon_f5.gtf")
    if [ "$mod" = noncode ]; then
        mkdir -p "$nadd"
        outs+=("$nadd/${prefix}_single_exon_f5.gtf" "$nadd/${prefix}_multi_exon_f5.gtf")
    fi
    for f in "${outs[@]}"; do : > "$f"; done
    awk -F'\t' -v OFS='\t' -v p="$outdir/$prefix" -v na="$nadd/$prefix" '
        FILENAME==ARGV[1]{
            if (FNR == 1) print "transcript_id\tchr\tstart\tend\tstrand\ttype\tn_exon\tlength\tFPKM\tsource\tstatus" > (p "_lncRNA_status.tsv")
            if ($11 == "") $11 = "pass"
            print > (p "_lncRNA_status.tsv")
            if ($10 == "known" && $11 != "known_lowFPKM") dest[$1] = p "_known_lncRNA_f1.gtf"
            if ($11 == "pass") {
                f = p "_" $6 "_exon_f5.gtf"
                if ($10 == "known") { dest2[$1] = na "_" $6 "_exon_f5.gtf" } else dest[$1] = f
            }
            next
        }
        /^#/{ next }
        {
            a = " " $9; k = index(a, " transcript_id \""); if (!k) next
            id = substr(a, k + 16); id = substr(id, 1, index(id, "\"") - 1)
            if (id in dest) print > dest[id]
            if (id in dest2) print > dest2[id]
        }' "$tmp/tx.tsv" "$input"

    local known="$outdir/${prefix}_known_lncRNA_f1.gtf" title="known"
    if [ "$mod" = noncode ]; then
        cat "$nadd/${prefix}_single_exon_f5.gtf" "$nadd/${prefix}_multi_exon_f5.gtf" > "$nadd/${prefix}_final_known.gtf"
        known="$nadd/${prefix}_final_known.gtf"; title="known NONCODE"
    fi
    cat "$known" "$outdir/${prefix}_single_exon_f5.gtf" "$outdir/${prefix}_multi_exon_f5.gtf" > "$outdir/${prefix}_final.gtf"

    awk -F'\t' -v t="$title" 'NR>1 && $11=="pass"{ if ($10=="known") k++; else n[$6]++ }
        NR>1 && $11=="known_pass"{ k++ }
        END{
            printf "\nTotal %d final lncRNA, out of these:\n", k + n["single"] + n["multi"]
            printf "%d are final single-exon lncRNA\n%d are final multi-exon lncRNA\n%d are %s lncRNA\n", n["single"], n["multi"], k, t
        }' "$outdir/${prefix}_lncRNA_status.tsv"
}

main(){
    echo "lncRNA.sh $VERSION: $input ($mod mode, $threads threads)"
    parse_gtf
    overlap_filters
    coding_filters
    write_outputs
    echo ""
    echo "Run succeed"
}

main
