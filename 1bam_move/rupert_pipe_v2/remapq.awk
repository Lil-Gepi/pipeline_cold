#!/Users/ychen/Library/bin/awk -f
# reestimate mapping qualities based on AS/XS

BEGIN{
	version="0.1"
	script_name="remapq.awk"
	pg_hdr_id="remapq"
	# smoothstep order to use
	if (NORD=="")
		NORD=7.5
	# maximum mapping quality to assign
	QMAX=60
	# store original mapping quality in this tag
	orig_mq_tag="XQ"
	# corresponding probability p = 10**(-q/10)
	p_qmax=10^(-QMAX/10)
	OFS="\t"
	write_pg_hdr=1
	pp_pg_id=""
	in_hdr_block=1
}
# find the last entry in the @PG block
/^@PG/{
	if (in_hdr_block) {
		for (i=2;i<=NF;i++) {
			split($i,x,":")
			if (x[1]=="ID") {
				pp_pg_id=x[2]
				break
			}
		}
	}
}
# output header lines unchanged
/^@/{
	if (in_hdr_block) {
		print
		next
	}
}
# append @PG line to header
write_pg_hdr{
	pgline=sprintf("@PG\tID:%s\tPN:%s\tVN:%s",pg_hdr_id,script_name,version)
	if (pp_pg_id!="")
		print pgline, "PP:" pp_pg_id
	else
		print pgline
	write_pg_hdr=0
	in_hdr_block=0
}
# process read lines (MUST come in pairs; TODO check header?)
{
	delete QNAME
	delete MAPQ
	delete NEWQ
	delete LINE
	for(i=1;i<=2;i++) {
		# read name
		QNAME[i]=$1
		# sanity check for pairs
		if (i==2 && QNAME[2]!=QNAME[1]) {
			print "Error in lines " NR-1 "/" NR ": file must be sorted by qname and only contain primary pairs!" > "/dev/stderr"
			exit
		}
		# original mapping quality
		MAPQ[i]=$5
		# find AS and XS fields in optional fields
		AS=0
		XS=0
		have_AS=0
		have_XS=0
		for (j=12;j<NF;j++) {
			split($j,x,":")
			switch (x[1]) {
				case "AS" :
					AS=int(x[3])
				       	have_AS=1
					break
				case "XS" :
					XS=int(x[3])
				       	have_XS=1
					break
			}
			if (have_AS && have_XS)
				break
		}
		if (XS>0)
			# heuristic mapping quality from AS/XS
			NEWQ[i]=p2q(1-smoothstepN(AS/(AS+XS)))
		else
			NEWQ[i]=QMAX
		# save line
		LINE[i]=$0
		if (i==1)
			# read mate line
			getline
	}
	# assign the same mapping quality to both members of the pair
	CONSQ=QMAX
	for (i in LINE) {
		if (MAPQ[i]<CONSQ)
			CONSQ=MAPQ[i]
		if (NEWQ[i]<CONSQ)
			CONSQ=NEWQ[i]
	}
	# modify output lines where necessary
	if (CONSQ<MAPQ[1] || CONSQ<MAPQ[2]) {
		 for (i in LINE) {
			 split(LINE[i],x)
			 # generate modified line
			 x[5]=CONSQ
			 rmvd_MQ=rmvd_PG=0
			 for (j=12;j<=NF;j++) {
				 split(x[j],y,":")
				 switch (y[1]) {
					 case "MQ" :
						 delete x[j]
						 rmvd_MQ=1
						 break
					 case "PG" :
						 delete x[j]
						 rmvd_PG=1
						 break
				 }
				 if (rmvd_MQ && rmvd_PG)
					 break
			 }
			 LINE[i]=sprintf("%s\tMQ:i:%i\t%s:i:%i\tPG:Z:%s",
			     collapse(x,"\t"),CONSQ,orig_mq_tag,MAPQ[i],pg_hdr_id)
		 }

	}
	print LINE[1]
	print LINE[2]
}
# rational smoothstep function from https://tpfto.wordpress.com/2019/03/28/on-a-rational-variant-of-smoothstep/
function smoothstepN(x) {
	# NORD is a global variable
	return 1/(1+((1-x)/x)^NORD)
	# this is the smoothstep polynomial of order 5
	#return (6*x^2-15*x+10)*x^3
}
# convert probability into phred quality
function p2q(p, q,r) {
	if (p>=1)
		q=0
	else if (p>p_qmax) {
		# q = -10*log10(p) = -10*log(x)/log(10) = -10*0.4342945*log(x)
		r=-4.342945*log(p)
		q=int(r)
		if (q<QMAX && r-q>=0.5) q++
	} else
		q=QMAX
	return q
}
function collapse(x,sep, i,r) {
	if (isarray(x)) {
		for (i in x)
			if (r)
				r=r sep x[i]
			else
				r=x[i]
	} else
		r=x
	return r
}
