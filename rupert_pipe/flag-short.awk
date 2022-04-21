#!/usr/local/bin/gawk -f
# flag templates shorter than the read length as QCFAIL

BEGIN{
	script_name="flag-short.awk"
	version="0.1"
	pg_hdr_id="flag-short"
	# check if READ_LENGTH was given on command line
	if (READ_LENGTH=="") {
		print "Warning: READ_LENGTH not set (use, e.g., -v READ_LENGTH=100), checking SEQ for each line ..." > "/dev/stderr"
		READ_LENGTH=0
		rl_from_seq=1
	} else {
		READ_LENGTH=int(READ_LENGTH)
		rl_from_seq=0
	}
	QCFAIL=512
	OFS="\t"
	write_pg_hdr=1
	pg_id_pp=""
	in_hdr_block=1
}
/^@PG/{
	if (in_hdr_block) {
		for (i=2;i<=NF;i++) {
			if ($i ~ /^ID:/) {
				split($i,x,":")
				pg_id_pp=x[2]
				break
			}
		}
	}
}
# output header unchanged
/^@/{
	print
	next
}
# write pg header
write_pg_hdr{
	PGLINE=sprintf("@PG\tID:%s\tPN:%s\tVN:%s\tCL:%s -v READ_LENGTH=%i",
	      pg_hdr_id,script_name,version,script_name,READ_LENGTH)
	if (pg_id_pp)
		print PGLINE, "PP:" pg_id_pp
	else
		print PGLINE
	write_pg_hdr=0
	in_hdr_block=0
}
# set READ_LENGTH from sequence length if requested
rl_from_seq{
	READ_LENGTH=length($10)
}
# process read lines
{
	TLEN=$9
	if (TLEN<0) TLEN*=-1
	if (TLEN>=READ_LENGTH) {
		print
		next
	}
	# adjust FLAG
	FLAG=or($2,QCFAIL)
	# generate & emit modified line
	LINE=$1 OFS FLAG
	for (i=3;i<=NF;i++) {
		if ($i ~ /^PG:/)
			continue
		else
			LINE=LINE OFS $i
	}
	print LINE, "PG:Z:" pg_hdr_id
}
