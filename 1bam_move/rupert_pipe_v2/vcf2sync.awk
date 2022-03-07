#!/usr/bin/gawk -f
# transform VCF into SYNC file
# # using
BEGIN{
	script_name="vcf2sync.awk"
	version="0.1"
	OFS="\t"
	# global state variables
	CHROM=""
	POS=""
	LINE=""
	# base indices in SYNC file: A:T:C:G:N:D
	SYX["A"]=1
	SYX["T"]=2
	SYX["C"]=3
	SYX["G"]=4
	SYX["N"]=5
	SYX["D"]=6
	# counters for deletions
	delete DEL
	delete DEL_REF
	delete DEL_RO
	DEL_CHROM=""
}
# ignore headers
/^#/{
	next
}
# process read lines
# we may encounter subsequent lines for the same position
$1==CHROM && $2==POS {
	warn_with_msg("duplicate position, merging with previous...")
	LINE=merge_lines(LINE,mk_sync_line())
	next
}
# emit stored DEL entries
DEL_CHROM {
	if (DEL_CHROM!=$1)
		# changed chromosome, emit all
		for (pos in DEL)
			emit_del_line(pos)
	else
		# same chrom, emit all with pos < current
		# should be ordered correctly as index is numeric
		for (pos in DEL)
			if (pos<$2)
				emit_del_line(pos)
	if (length(DEL)==0)
		DEL_CHROM=""
}
# emit stored line and make a new one
{
	if (LINE) {
		print LINE
		LINE=""
	}
	LINE=mk_sync_line()
}
# emit last line when we reach the end
END{
	if (LINE)
		print LINE
}
function mk_sync_line( chrom,pos,ref,alt,x,info,del_len,real_ref,i,format,adx,count,y,line) {
	# fields
	chrom=$1
	pos=$2
	ref=$4
	if (ref==".")
		skip_with_msg("REF information missing")
	alt=$5
	# only allow bi-allelic sites
	split(alt,x,",")
	if (length(x)>1)
		skip_with_msg("multi-allic sites not allowed")
	info=$8
	# only allow SNP sites, check for fake REF
	split(info,x,";")
	real_ref=ref
	del_len=0
	for (i in x)
	switch (x[i]) {
		case "INDEL" : 
			del_len=length(ref)-length(alt)
			if (del_len<0)
				skip_with_msg("ignoring INS")
			break
		case "RMOD"  :
			warn_with_msg("real REF unknown, setting to '.'")
			real_ref="."
			break
	}
	# handle SNP
	# double check that we have single bases
	if (del_len==0 && (length(ref)>1 || length(alt)>1))
		skip_with_msg("complex variants not supported")
	format=$9
	# position of AD
	split(format,x,":")
	adx=0
	for (i in x)
		if (x[i]=="AD") {
			adx=i
			break
		}
	if (adx==0)
		skip_with_msg("FORMAT/AD tag not found")
	# write CHROM POS REF [SAMPLE]
	line=chrom OFS pos OFS real_ref
	if (del_len==0) {
		# handle SNP
		for(i=10;i<=NF;i++) {
			for (j in SYX)
				count[SYX[j]]=0
			split($i,x,":")
			split(x[adx],y,",")
			count[SYX[ref]]=int(y[1])
			count[SYX[alt]]=int(y[2])
			if (DEL[pos][i]) {
				count[SYX["D"]]=DEL[pos][i]
				delete DEL[pos][i]
			}
			line=line OFS collapse(count,":")
		}
		delete DEL[pos]
		delete DEL_REF[pos]
		if (length(DEL)==0)
			DEL_CHROM=""
		CHROM=chrom
		POS=pos
		return line
	} else {
		# handle DEL
		split(ref,x,"")
		for (j=1;j<=del_len;j++) {
			DEL_REF[pos+j]=x[1+j]
		}
		for(i=10;i<=NF;i++) {
			split($i,x,":")
			if (x[adx]==".") {
				y[1]=0
				y[2]=0
			} else
				split(x[adx],y,",")
			for (j=1;j<=del_len;j++) {
				if (DEL[pos+j][i]) {
					DEL_RO[pos+j][i]+=int(y[1])
					DEL[pos+j][i]+=int(y[2])
				} else {
					DEL_RO[pos+j][i]=int(y[1])
					DEL[pos+j][i]=int(y[2])
				}
			}
		}
		DEL_CHROM=chrom
	}
}
function merge_lines(L1,L2, u,v,x,y,i,j,n) {
	split(L1,x)
	split(L2,y)
	# check sample number
	if (length(y)!=length(x))
		skip_with_msg("duplicate position, but different sample number?")
	# check REF
	if (x[3]!=y[3])
		skip_with_msg("REFs do not match") 
	n=length(x)
	for (i=4;i<=n;i++) {
		split(x[i],u,":")
		split(y[i],v,":")
		for (j in u)
			u[j]+=v[j]
		x[i]=collapse(u,":")
	}
	return collapse(x,"\t")
}
function skip_with_msg(m) {
	print "Skipping record " NR ": " m > "/dev/stderr"
	next
}
function warn_with_msg(m) {
	print "Issue in record " NR ": " m > "/dev/stderr"
}
function collapse(x,sep, i,r) {
	if (isarray(x)) {
		for (i in x)
			if (r)
				r=r sep x[i]
			else
				r=x[i] ""
	} else
		r=x
	return r
}
function emit_del_line(pos, line,i,count) {
	line=DEL_CHROM OFS pos OFS DEL_REF[pos]
	for (i in DEL[pos]) {
		for (j in SYX)
			count[SYX[j]]=0
		count[SYX[DEL_REF[pos]]]=DEL_RO[pos][i]
		count[SYX["D"]]=DEL[pos][i]
		line=line OFS collapse(count,":")
	}
	print line
	delete DEL[pos]
	delete DEL_REF[pos]
	delete DEL_RO[pos]
}
