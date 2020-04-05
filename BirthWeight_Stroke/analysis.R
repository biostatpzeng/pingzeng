#LDSC
	bash
	for dx in EGG_BirthWeight_Fetal_Effect_2019 EGG_BirthWeight_Maternal_Effect_2019; do
	for ex in AIS CES LAS SVS; do
	ldsc=/public/software/apps/ldsc/ldsc.py
	export PATH=/public/software/apps/anaconda2/bin:$PATH
	python ${ldsc} \
	--rg ${dx}.txt,${ex}.txt \
	--ref-ld-chr ${home}/EUR_LDscore/ \
	--w-ld-chr   ${home}/EUR_LDscore/ \
	--out ./${ex}_${dx}
	done
	done


#	Summary of Genetic Correlation Results
#	p1                                   p2      rg      se       z       p  h2_obs  h2_obs_se  h2_int  h2_int_se  gcov_int  gcov_int_se
#	EGG_BirthWeight_Maternal_Effect_2019.txt  IschemicStroke_AS_EUR_Malik2018.txt -0.1763  0.0539 -3.2723  0.0011  0.0107     0.0009  1.0153     0.0044   -0.0003       0.0029



# MR
	library(data.table)
	library(MendelianRandomization)
	bx = data.frame(fread(paste("index_snp_EGG_BirthWeight_Fetal_2019_5e8.csv",sep=""),head=T))
	bx = data.frame(fread(paste("index_snp_EGG_BirthWeight_maternal_2019_5e8.csv",sep=""),head=T))
	bx=na.omit(bx)
	x = c("AIS","CES","LAS","SVS")
	m = length(x)
	rex = matrix(NA, m, 6)
	snp = as.character(bx$SNP)
	for (j in 1:length(x)){
		try({
			try({
				mx0 = data.frame(fread(paste(x[j],".txt",sep=""),sep='\t',head=T))
				mx = mx0[which(as.character(mx0$SNP)%in%snp),]
			})
			try({
				yx = merge(bx,mx,by.x="SNP",by.y="SNP",all.x=T)
				index1=which(toupper(yx$INC_ALLELE.x)==toupper(yx$INC_ALLELE.y))
				if (length(index1)<length(yx$INC_ALLELE.x)){yx$BETA.y[-index1]=yx$BETA.y[-index1]*(-1)}
				yx=na.omit(cbind(
				as.numeric(paste(yx$BETA.x)),
				as.numeric(paste(yx$SE.x)),
				as.numeric(paste(yx$BETA.y)),
				as.numeric(paste(yx$SE.y)),
				as.numeric(paste(yx$P.y))
				))
				colnames(yx)=c("BETA.x","SE.x","BETA.y","SE.y","P.y")
				yx=as.data.frame(yx)
				index=which(yx$SE.y>0)
				yx=as.data.frame(yx[index,])
				index2=which(yx$P.y<=(0.05/length(yx$P.y)))
				if (length(index2)>0) {yx=yx[-index2,]}
				fit = mr_ivw(mr_input(bx = yx$BETA.x, bxse = yx$SE.x, by = yx$BETA.y, byse = yx$SE.y),model="fixed")
				fit = mr_ivw(mr_input(bx = yx$BETA.x, bxse = yx$SE.x, by = yx$BETA.y, byse = yx$SE.y),model="random")
				fit = mr_median(mr_input(bx = yx$BETA.x, bxse = yx$SE.x, by = yx$BETA.y, byse = yx$SE.y))
				fit = mr_egger(mr_input(bx = yx$BETA.x, bxse = yx$SE.x, by = yx$BETA.y, byse = yx$SE.y))
			})
			rm(fit)
			rm(yx)
			print(j)
		})
	}
	cc = c("N","beta_fix","beta_se_fix","p_value_fix","Heter.Stat","Heter.Stat.p","beta_random","beta_se_random","p_value_random","beta_median","beta_se_median","p_value_median","beta_egger","beta_se_egger","p_value_egger","Intercept_egger","Intercept_se_egger","Intercept.p_egger")
	colnames(rex) = cc
	rownames(rex) = x

#	> rex_maternal
#					   N   beta_fix beta_se_fix  p_value_fix
#	IschemicStroke_AS_EUR_Malik2018  238 -0.1166151  0.04468245 0.0090578521
#	IschemicStroke_CES_EUR_Malik2018 241 -0.1449659  0.09240212 0.1166808010
#	IschemicStroke_LAS_EUR_Malik2018 240 -0.4041020  0.12022815 0.0007762544
#	IschemicStroke_SVS_EUR_Malik2018 239 -0.3823822  0.11229824 0.0006614951
#					 Heter.Stat Heter.Stat.p beta_random
#	IschemicStroke_AS_EUR_Malik2018    381.5418 7.619406e-09  -0.1166151
#	IschemicStroke_CES_EUR_Malik2018   331.4476 8.261889e-05  -0.1449659
#	IschemicStroke_LAS_EUR_Malik2018   316.4655 5.782664e-04  -0.4041020
#	IschemicStroke_SVS_EUR_Malik2018   286.0466 1.788967e-02  -0.3823822
#					 beta_se_random p_value_random beta_median
#	IschemicStroke_AS_EUR_Malik2018      0.05669358    0.039692365  -0.1289862
#	IschemicStroke_CES_EUR_Malik2018     0.10858848    0.181875526  -0.1087780
#	IschemicStroke_LAS_EUR_Malik2018     0.13834721    0.003489928  -0.4951629
#	IschemicStroke_SVS_EUR_Malik2018     0.12311271    0.001896671  -0.2921061
#					 beta_se_median p_value_median beta_egger
#	IschemicStroke_AS_EUR_Malik2018       0.0728579     0.07666372 -0.1402764
#	IschemicStroke_CES_EUR_Malik2018      0.1449469     0.45297291 -0.1785304
#	IschemicStroke_LAS_EUR_Malik2018      0.1971142     0.01200283 -0.5694126
#	IschemicStroke_SVS_EUR_Malik2018      0.1680629     0.08219738 -0.5902531
#					 beta_se_egger p_value_egger Intercept_egger
#	IschemicStroke_AS_EUR_Malik2018     0.08361112   0.093400994    0.0005082196
#	IschemicStroke_CES_EUR_Malik2018    0.16229545   0.271317499    0.0007232092
#	IschemicStroke_LAS_EUR_Malik2018    0.20626035   0.005768648    0.0035445238
#	IschemicStroke_SVS_EUR_Malik2018    0.18274409   0.001238147    0.0044410207
#					 Intercept_se_egger Intercept.p_egger
#	IschemicStroke_AS_EUR_Malik2018         0.001350091         0.7065946
#	IschemicStroke_CES_EUR_Malik2018        0.002630006         0.7833286
#	IschemicStroke_LAS_EUR_Malik2018        0.003319466         0.2856111
#	IschemicStroke_SVS_EUR_Malik2018        0.002930057         0.1296010
#
#	> rex_fetal
#					   N    beta_fix beta_se_fix p_value_fix
#	IschemicStroke_AS_EUR_Malik2018  239 -0.03967466  0.03796759  0.29604088
#	IschemicStroke_CES_EUR_Malik2018 243  0.01027267  0.07886842  0.89636809
#	IschemicStroke_LAS_EUR_Malik2018 242 -0.16731417  0.10124064  0.09840447
#	IschemicStroke_SVS_EUR_Malik2018 241 -0.16442884  0.09448895  0.08182544
#					 Heter.Stat Heter.Stat.p beta_random
#	IschemicStroke_AS_EUR_Malik2018    380.1659 1.272400e-08 -0.03967466
#	IschemicStroke_CES_EUR_Malik2018   324.5798 3.072518e-04  0.01027267
#	IschemicStroke_LAS_EUR_Malik2018   324.7407 2.559382e-04 -0.16731417
#	IschemicStroke_SVS_EUR_Malik2018   303.9354 3.234934e-03 -0.16442884
#					 beta_se_random p_value_random beta_median
#	IschemicStroke_AS_EUR_Malik2018      0.04798562      0.4083487 -0.06329793
#	IschemicStroke_CES_EUR_Malik2018     0.09133899      0.9104527  0.09515173
#	IschemicStroke_LAS_EUR_Malik2018     0.11752079      0.1545338 -0.14204247
#	IschemicStroke_SVS_EUR_Malik2018     0.10633249      0.1220165 -0.17539463
#					 beta_se_median p_value_median  beta_egger
#	IschemicStroke_AS_EUR_Malik2018       0.0646122      0.3272544 -0.03297299
#	IschemicStroke_CES_EUR_Malik2018      0.1308420      0.4670873  0.01943201
#	IschemicStroke_LAS_EUR_Malik2018      0.1659774      0.3921118 -0.28526867
#	IschemicStroke_SVS_EUR_Malik2018      0.1595905      0.2717552 -0.15623446
#					 beta_se_egger p_value_egger Intercept_egger
#	IschemicStroke_AS_EUR_Malik2018     0.09256952     0.7216930   -0.0001479578
#	IschemicStroke_CES_EUR_Malik2018    0.17520270     0.9116864   -0.0002021162
#	IschemicStroke_LAS_EUR_Malik2018    0.22593302     0.2067240    0.0026155874
#	IschemicStroke_SVS_EUR_Malik2018    0.20656149     0.4494345   -0.0001803570
#					 Intercept_se_egger Intercept.p_egger
#	IschemicStroke_AS_EUR_Malik2018         0.001748375         0.9325588
#	IschemicStroke_CES_EUR_Malik2018        0.003297931         0.9511316
#	IschemicStroke_LAS_EUR_Malik2018        0.004277222         0.5408584
#	IschemicStroke_SVS_EUR_Malik2018        0.003897460         0.9630906

# Multivariable MR
	library(data.table)
	library(MendelianRandomization)
	x=c("AIS","CES","LAS","SVS")
	ex=c("Early_GWG","Late_GWG","Total_GWG","EGG_Fetal_gest_duration_NComms2019","EGG_BMIchildhood","EGG_growth_PG",
		"EGG_TannerStage","EGG_HeadCircumference","EGG_CHILDHOOD_OBESITY_EUR_2019",
		"EGG_growth_PT","EGG_growth_1012","EGG_BirthLength","SSGAC_AgeFirstBirth",
		"Alcohol_continuous","CCACE_Income_UKBHill2016","BMI_Locke_UKBiobank_2018",
		"Bodyfat","CAD_add_2015","College","FastingGlucose","FastingInsulin",
		"FHS_weight","HbA1c_METAL_European","HDL2013","LDL2013","TC2013",
		"TG2013","HeartRate","Height_Yang2012","HEIGHT_Wood","HEIGHT_Randall2013",
		"HIP_COMBINED_EUR","HIPadjBMI_COMBINED_EUR","HOMAB","HOMAIR",
		"hrGlucose_AdjustedForBMI","Obesity1","Overweight",
		"SSGAC_EducationalAttainment_PNAS2014","T2D_BMI_2017","W1_hypertension",
		"WCadjBMI_COMBINED_EUR","WC_Randall2013","WEIGHT_Randall2013","WHRadjBMI_Randall2013",
		"WHR_Randall2013","GSCAN_CigarettesPerDay","GSCAN_DrinksPerWeek",
		"GSCAN_SmokingCessation","GSCAN_AgeOfInitiation","GSCAN_SmokingInitiation")
	m = length(ex)
	rex = matrix(NA, m, 5)
	for (i in 1:m){
		yx0 = data.frame(fread(paste("AIS_multi_M.txt",sep=""),sep='\t',head=T))
		snp=yx0$SNP
		x0 = data.frame(fread(paste(ex[i],".txt",sep=""),sep='\t',head=T))
		ax1=x0[which(x0$SNP%in%snp),]
		yx1=merge(yx0,ax1,by="SNP",all.x=T)
		index1=which(toupper(yx1$INC_ALLELE.x)==toupper(yx1$INC_ALLELE))
		if (length(index1)<length(yx1$INC_ALLELE.x)){yx1$BETA[-index1]=yx1$BETA[-index1]*(-1)}
		yx=na.omit(yx1)
		yy=lm(yx$BETA.y~yx$BETA-1,weights=yx$SE.y^2)$residuals
		fit=lm(yy~yx$BETA.x-1,weights=yx$SE.y^2)
		rex[i,1:5] = c(length(yy),summary(fit)$coef)
		print(i)
	}
#	variables	N	beta	se	t	p
#	EGG_Fetal_gest_duration_NComms2019	236	-0.237627913045698	0.0755867818409527	-3.14377603144565	0.00188257884900272
