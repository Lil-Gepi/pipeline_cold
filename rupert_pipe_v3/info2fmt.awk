#!/usr/local/bin/gawk -f
#
# copy INFO tags (specifically DP and MQ0) to the genotype field
#
# Rupert Mazzucco, 2022

BEGIN{
	version="0.0.9001"
	script_name="info2fmt"
	OFS="\t"
	# index in gentoype field
	delete FMT
	delete GTX
	# info TAGS to copy over
	delete TAGS
	# can be passed in via the "tags" variable
	if (tags) {
		split(tags,x,",")
		for (i in x)
			TAGS[x[i]]=1
	} else {
		TAGS["DP"]=1
		TAGS["MQ0"]=1
	}
	# header entries to add
	delete INFOHDR
	delete FMTHDR
	# make sure we traverse arrays in order after split
	PROCINFO["sorted_in"]="@ind_num_asc"
	write_hdr_line=1
}
# build FORMAT headers corresponding to INFO headers
/^##INFO/{
	split($0,x,"(=|,)")
	id=x[3]
	INFOHDR[id]=$0
	if (id in TAGS) {
		line=$0
		sub("^##INFO","##FORMAT",line)
		FMTHDR[id]=line
	}
	next
}
# delete existing FORMAT header where we have a new one
/^##FORMAT/{
	split($0,x,"(=|,)")
	id=x[3]
	if (!(id in TAGS))
		FMTHDR[id]=$0
	next
}
# emit other headers unchanged
/^##/{

	print
	next
}
# emit our headers
write_hdr_line{
	for (i in TAGS) {
		if (!(i in INFOHDR))
			print "Error: no header found for tag INFO/" i ", exiting." > "/dev/stderr"
		if (!(i in FMTHDR))
			print "Error: failed to construct header for  FORMAT/" i ", exiting." > "/dev/stderr"
	}
	PROCINFO["sorted_in"]="@ind_str_asc"
	#asort(INFOHDR)
	for (i in INFOHDR)
		print INFOHDR[i]
	#asort(FMTHDR)
	for (i in FMTHDR)
		print FMTHDR[i]
	PROCINFO["sorted_in"]="@ind_num_asc"
	print "##poptools_Command=" script_name "; Version=" version "; CommandLine=" command_line "; Date=" strftime()
	write_hdr_line=0
}
/^#CHROM/{
	print
	next
}
# basic VCF fields (ignore ID, QUAL, and FILTER)
{
	CHROM=$1
	POS=$2
	REF=$4
	ALT=$5
	INFO=$8
	# ids and vals of info tags to copy to genotype field
	delete i2f
	split(INFO,x,";")
	for (i in x) {
		split(x[i],y,"=")
		if (y[1] in TAGS)
			i2f[y[1]]=y[2]
	}
}
# rebuild the FORMAT index if FORMAT changes (unlikely, but possible)
$9!=FORMAT{
	FORMAT=$9
	split(FORMAT,FMT,":")
	delete GTX
	for (i in FMT)
		GTX[FMT[i]]=i
}
# build new format & genotype fields and print
{
	delete fmt
	delete gt
	for (i in FMT)
		fmt[i]=FMT[i]
	split($10,gt,":")
	for (i in i2f) {
		if (i in GTX) {
			gt[GTX[i]]=i2f[i]
		} else {
			fmt[length(fmt)+1]=i
			gt[length(gt)+1]=i2f[i]
		}
	}
	print CHROM, POS, ID, REF, ALT, ".", ".", INFO, collapse(fmt,":"), collapse(gt,":")
}
# collapse an ordered array with consecutive indices 1..length(x) to a string
function collapse(x,sep, i,S) {
	if (isarray(x)) {
		S=x[1]
		for (i=2;i<=length(x);i++)
			S=S sep x[i]
	} else
		S=x
	return S
}
