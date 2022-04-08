#!/usr/local/bin/gawk -f

# TODO
# 1.h annotate QUAL based on mapping qualities, base qualities
# 1.h.1 use coverage/quality anomalies (based on INFO/DP-sum(FMT/DP))
# 1.h.2 use quality anomalies (based on abs(INFO/QR / INFO/RO - INFO/QA / INFO/AO))
# 1.i make GT call to heuristically optimize QUAL (?)

BEGIN{
	# script version
	version="0.1"
	script_name="post-merging.awk"

	OFS="\t"
	in_format=0
	in_info=0
	after_hdr=0
	fmt_done=0
	inf_done=0
}
/^##ALT=<ID=\*/{
	# don't need <*> allele anymore
	next
}
/^##INFO/{
	in_info=1
	print
	next
}
/^##FORMAT/{
	in_format=1
	print
	next
}
/^#/{
	if (inf_done==0 && in_info==1) {
		in_info=0
		print "##INFO=<ID=N_ALT,Number=1,Type=Integer,Description=\"Number of ALT alleles\">"
		print "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Raw read depth\">"
		print "##INFO=<ID=RO,Number=1,Type=Integer,Description=\"Number of high-quality REF counts\">"
		print "##INFO=<ID=AO,Number=A,Type=Integer,Description=\"Number of high-quality ALT counts\">"
		print "##INFO=<ID=QR,Number=1,Type=Integer,Description=\"Sum of phred-scaled REF qualities\">"
		print "##INFO=<ID=QA,Number=A,Type=Integer,Description=\"Sum of phred-scaled ALT qualities\">"
		print "##INFO=<ID=RMOD,Number=0,Type=Flag,Description=\"REF not observed and replaced by ALT{1}\">"
		inf_done=1
	}
	if (fmt_done==0 && in_format==1) {
		in_format=0
		print "##FORMAT=<ID=SAC,Number=1,Type=Integer,Description=\"Sum of allele counts, SAC = sum(AD)\">"
		print "##FORMAT=<ID=AF,Number=A,Type=Float,Description=\"Observed ALT frequencies\">"
		print "##FORMAT=<ID=XF,Number=A,Type=Float,Description=\"Estimated ALT frequencies\">"
		fmt_done=1
	}
	if ($0 ~ /^#CHROM/)
		print "##poptools_annotateCommand=post-merging.awk; Version=" version "; Date=" strftime()
	print
	next
}
{
	CHROM=$1
	POS=$2
	ID=$3
	REF=$4
	ALT=$5
        # 1.a drop ALT=<*> positions
	if (ALT=="<*>")
		next
        # and alleles (also fix FMT/AD and FMT/QS)
	split(ALT,x,",")
	j=0
	fix_adqs=0
	to_del=0
       	for (i in x)
	       	if (x[i]=="<*>") {
			fix_adqs=1
			delete x[i]
			# ALT has length A, but AD+QS have R, i.e., A+1
			to_del = i+1
			# to_del now contains the R-index to remove
			break
		}
        # 1.b annotate positions with INFO/N_ALT
	N_ALT=length(x)
	ALT=collapse(x,",")
	QUAL=$6
	FILTER=$7
	INFO=$8
	split(INFO,x,";")
	TYPE="snp"
	for (i in x)
		if (x[i]=="INDEL")
			TYPE="indel"
	FORMAT=$9
	# and RO,QR,AO,QA from sample AD and QS
	# indices of FMT/AD and FMT/QS
	adx=0
	qsx=0
	split(FORMAT,x,":")
	for (i in x)
		switch (x[i]) {
			case "AD" :
				adx=i
				break
			case "QS" :
				qsx=i
				break
			case "DP" :
				dpx=i
				break
		}
	# will add AF and XF annotations
	sacx=length(x)+1
	afx=length(x)+2
	xfx=length(x)+3
	x[sacx]="SAC"
        x[afx]="AF"
	x[xfx]="XF"
	FORMAT=collapse(x,":")
	# go over samples and build up various INFO tags
	# fix up AD and QS where necessary
	# add AF and XF
	DP=0
	RO=0
	QR=0
	delete AO
	delete QA
	for (i=1;i<=N_ALT;i++) {
		AO[i]=0
		QA[i]=0
	}
	for (i=10;i<=NF;i++) {
		split($i,x,":")
		if (x[dpx]==".")
			x[dpx]=0
		DP+=x[dpx]
		split(x[adx],ad,",")
		split(x[qsx],qs,",")
		if (length(ad)==1) {
			# samples without info have only "." instead of an R-array ...
			# recreate the R-array, taking into account whether or not
			# one will be trimmed off. QS will generally conform with AD.
			if (fix_adqs==1)
				dn=2
			else
				dn=1
			for (j=1;j<=N_ALT+dn;j++) {
				ad[j]="."
				qs[j]="."
			}
		}
		for (j in ad) {
			if (ad[j]==".") ad[j]=0
			if (qs[j]==".") qs[j]=0
		}
		if (fix_adqs==1) {
			delete ad[to_del]
			delete qs[to_del]
		}
		RO+=ad[1]
		adsum=ad[1]
		QR+=qs[1]
		k=1
		for (j in ad) {
			if (j==1) continue
			adsum+=ad[j]
			AO[k]+=ad[j]
			QA[k]+=qs[j]
			k++
		}
		x[adx]=collapse(ad,",")
		x[qsx]=collapse(qs,",")
		# observed and estimated allele fractions, AF and XF
		delete af
		delete xf
		xfsum=adsum+length(ad)
		k=1
		for (j in ad) {
			if (j==1) continue
			if (adsum>0)
				af[k]=ad[j]/adsum
			else
				af[k]="."
			xf[k]=(ad[j]+1)/xfsum
			k++
		}
		# add to genotype vector
		x[sacx]=adsum
		x[afx]=collapse(af,",")
		x[xfx]=collapse(xf,",")
		# recompose genotype field
		$i=collapse(x,":")
	}
	# 1.c flag RO=0 positions and change REF to ALT{1}, ALT{1} to ALT{2} etc.
	mod_ref=0
	if (RO==0 && N_ALT>1 && TYPE=="snp") {
		split(ALT,x,",")
		REF=x[1]
		delete x[1]
		ALT=collapse(x,",")
		# 1.d adjust all R and A tags, specifically INFO/{RO,QR,AO,QA}
		RO=AO[1]
		delete AO[1]
		QR=QA[1]
		delete QA[1]
		# and N_ALT
		N_ALT--;
		# and FMT/{AD,QS,AF,XF}
		for (i=10;i<=NF;i++) {
			split($i,x,":")
			split(x[adx],ad,",")
			delete ad[1]
			x[adx]=collapse(ad,",")
			split(x[qsx],qs,",")
			delete qs[1]
			x[qsx]=collapse(qs,",")
			split(x[afx],af,",")
			delete af[1]
			x[afx]=collapse(af,",")
			# recalculate XF because number of classes is now -1
			adsum=0
			for (j in ad)
				adsum+=ad[j]
			xfsum=adsum+length(ad)
			delete xf
			for (j=1;j<=N_ALT;j++)
				xf[j]=(ad[j+2]+1)/xfsum
			x[xfx]=collapse(xf,",")
			$i=collapse(x,":")
		}
		# flag modified REF
		mod_ref=1
	}
	# TODO make call, assess QUAL ...
	QUAL="."
	FILTER="."
	# update INFO field
	INFO_NEW=sprintf("N_ALT=%i;DP=%i;RO=%i;AO=%s;QR=%i;QA=%s",
		N_ALT,DP,RO,collapse(AO,","),QR,collapse(QA,","))
	if (INFO==".")
		INFO=INFO_NEW
	else
		INFO=INFO_NEW ";" INFO
	if (mod_ref==1)
		INFO=INFO ";RMOD"
	# emit line
	LINE=sprintf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s",
	    CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT)
	for (i=10;i<=NF;i++)
		LINE=sprintf("%s\t%s",LINE,$i)
	print LINE
}
function collapse(x,sep, i,r) {
	if (isarray(x)) {
		for (i in x)
			if (r) r=r sep x[i]
		        else r=x[i] ""
		return r
	} else
		return x
}
