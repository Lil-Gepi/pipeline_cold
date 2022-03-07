#!/bin/bash
# Rupert Mazzucco, 2021

version="0.1"

usage="
  $0 in.bcf [out.bcf]
  $0 -h
  $0 -v
  
  Preprocess a single-sample VCF (in any format supported by bcftools) as produced
  by bcftools mpileup by turning most INFO tags into FORMAT tags. Default output is 
  to stdout in VCF format. Specify '-' to get uncompressed BCF output to stdout.
"

if [[ -z "$*" || "$1" = "-h" ]]; then
	echo "$usage"
	exit
fi

if [[ "$1" = "-v" ]]; then
	echo "$version"
	exit
fi

ivcf="$1"
ovcf=""

if [[ -n "$2" ]]; then
	case "$2" in
		*.vcf.gz) ovcf="-Oz -o $2" ;;
		*.vcf)    ovcf="-Ov -o $2" ;;
	        *.bcf)    ovcf="-Ob -o $2" ;;
		-)        ovcf="-Ou -o -" ;;
		*)        echo "Unknown output format requested, exiting ..." && exit 1 ;;
	esac
fi

trap 'rm -rf $tdir' EXIT

# temporary files
TMPDIR=${TMPDIR:-/tmp}
tdir=$(mktemp -t -d precall.XXXXX)
ann="$tdir/ann.txt.gz"
hdr="$tdir/ann.hdr"

# INFO tags to extract and transform into FORMAT tags
itags="
DP	1	Integer	Raw read depth
RPBZ	1	Float	Mann-Whitney U-z test of Read Position Bias (closer to 0 is better)
MQBZ	1	Float	Mann-Whitney U-z test of Mapping Quality Bias (closer to 0 is better)
BQBZ	1	Float	Mann-Whitney U-z test of Base Quality Bias (closer to 0 is better)
MQSBZ	1	Float	Mann-Whitney U-z test of Mapping Quality vs Strand Bias (closer to 0 is better)
SCBZ	1	Float	Mann-Whitney U-z test of Soft-Clip Length Bias (closer to 0 is better)
FS	1	Float	Phred-scaled p-value using Fisher's exact test to detect strand bias
MQ0F	1	Float	Fraction of MQ0 reads (smaller is better)
I16	16	Float	Auxiliary tag, see bcftools source: bcf_callret1_t in bam2bcf.h
"

echo "$itags" | while read id num type desc; do
         test -z "$id" && continue
         printf "##FORMAT=<ID=%s,Number=%s,Type=%s,Description=\"%s\">\n" "$id" "$num" "$type" "$desc"
done > "$hdr"

# fields and INFO tags to extract
rmtags=""
addtgs="CHROM,POS,REF,ALT,"
printf -v fmtstr '%%%s\\t' CHROM POS REF ALT
while read tag rest; do
         test -z "$tag" && continue
	 fmtstr="${fmtstr}%$tag\\t"
	 rmtags="${rmtags}INFO/$tag,"
	 addtgs="${addtgs}FORMAT/$tag,"
done < <(echo "$itags")
fmtstr=${fmtstr%t}n
rmtags="INFO/VDB,INFO/SGB,INFO/QS,${rmtags}FORMAT/PL,FORMAT/DP"
addtgs=${addtgs%,}

CMD="(bcftools query -f \"$fmtstr\" \"$ivcf\" | bgzip > \"$ann\") && tabix -s1 -b2 -e2 \"$ann\""
eval "$CMD"

CMD="bcftools annotate -x $rmtags -Ou \"$ivcf\" | bcftools annotate -a \"$ann\" -c $addtgs -h \"$hdr\" $ovcf -"
eval "$CMD"
