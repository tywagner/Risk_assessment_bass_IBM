options(scipen=999)

###########################simulation parameter
###########################baseline setting: temp1, flow1, chem.on=0 (no chem effect)
replicate=10
sim.yr=200
river=10*2000*5 #wide*long*deep (m3)
n.initial=1000 #initial size
stru.initial=read.table("input_ini_size_structure.txt",header=T, sep='\t',row.names=NULL) #read in iniital size structure;1st col-lower bound,2nd-upper,3rd-proportion
sex.female=0.5 #sex ratio 1:1

############################ biology module
########growth
swimup.size.low=8.5 #mm; lower bound of size at swimming up
swimup.size.up=9.5 #mm; upper bound of size at swimming up
disper.size=20 #mm; size at dispersal
juv.size.low=90 #mm size becoming juvenile
max.size=560 #mm
linf=797 #mm
k.mean.normal=0.008 #per month; k mean under normal condition
logk.cv=0.2
logk.sd=round(abs(logk.cv*log(k.mean.normal)),2) #k sd at log-scal, 0.97

dt=1 #time step 1 month

#########survival
mc=0.034 #constant natural mortality, per month; Then et al. 2015 based on max age of 15 in Susquehanna River
m0=round((max.size-swimup.size.low)*mc/(log(max.size)-log(swimup.size.low)),3)
# egg.viable.prop=0.34 #viable egg proportion, i.e., survived egg propotion among deposited/produced egg; Raffetto et al 1990
# egg.fry.surv=0.22 #proportion of viable egg surving to larvae/fry; Raffetto et al. 1990
egg.surv.normal=0.002 #=egg.viable.prop*egg.fry.surv

############reproduction
# sp.season=4:6 #spawning season Apr-June
mat.size.low=175 #maturity: fish<175mm 0, ≥225mm 1, between=linear increase; Hoopes 1987
mat.size.up=225

nest.surv.normal=0.574 #survival of nest related to parental care of guarding male and environmental condition
fert=0.7 #fertilization sucess
log.ddays.int=17.669 #intercept for calculating log.ddays(cumulative degree days before a male starts nesting);DeAngelis et al 1991
log.ddays.slope=-3.484 #slope for calculating log.ddays(cumulative degree days before a male starts nesting)
fec.int=-12251.5 #intercept for calculating fecundity
fec.slope=59.4 #slope for calculating fecundity

##################################effect module
###########density dependence
max.density=0.8 #ind/m3
#scalar for egg survival: 1 before 1/2 max.density, linearly decreasing to 0 at max.density
#scalar for fry mortality: 1 before 1/2 max.density, lnearly increasing to 2 at max.density

############temperature effect
temp.data=read.table("input_temp.txt",header=T, sep='\t',row.names=NULL) #monthly temp time series
temp=temp.data$temp1 #switch between temp scenarios

ddays.low.data=read.table("input_ddays_low.txt",header=T, sep='\t',row.names=NULL) #lower bound of monthly ddays
ddays.low=ddays.low.data$temp1 #switch between temp scenarios

ddays.up.data=read.table("input_ddays_up.txt",header=T, sep='\t',row.names=NULL) #upper bound of monthly ddays
ddays.up=ddays.up.data$temp1 #switch between temp scenarios

nest.surv.temp=c(0,13,21,38) #ºC, spawning cease when temp<13 or >21 ºC
growth.temp=c(0,25,29,35) #ºC, growth cease when freezing or >35 ºC,optimal 25-29
egg.surv.temp=c(0,13,25,38) #ºC, survival cease when freezing or >38 ºC,optimal 13-25
fry.surv.temp=c(0,20,29,38) #ºC, survival cease when freezing or >38 ºC,optimal 20-29
juv.surv.temp=c(0,28,31,35) #ºC, survival cease when freezing or >35 ºC,optimal 28-31
adu.surv.temp=c(0,21,27,32) #ºC, survival cease when freezing or >32 ºC,optimal 21-27

#function to calculate temp scalar for mortality of fry, juvenile and adult given temp and temp threasholds;scalar=1 for optimal temp, linearly increase to 2 at upper and lower lethal temps
fun.temp.scalar.mort=function(tt.vec,tt){
	if(tt>=tt.vec[2] & tt<=tt.vec[3]){
		scalar=1
	}
	if(tt>tt.vec[1] & tt<tt.vec[2]){
		scalar=1+((2-1)*(tt.vec[2]-tt)/(tt.vec[2]-tt.vec[1]))
	}
	if(tt>tt.vec[3] & tt<tt.vec[4]){
		scalar=1+((2-1)*(tt-tt.vec[3])/(tt.vec[4]-tt.vec[3]))
	}
	if(tt<=tt.vec[1] | tt>=tt.vec[4]){
		scalar=10 #assign a value of 10 is to ensure survival=0 when temp beyond the survival range
	}
	scalar
}

#function to calculate temp scalar for survival of nest and eggs given temp and temp threasholds;scalar=1 for optimal temp, linearly decrease to 0 at upper and lower lethal temps
fun.temp.scalar.surv=function(tt.vec,tt){
	if(tt>=tt.vec[2] & tt<=tt.vec[3]){
		scalar=1
	}
	if(tt>tt.vec[1] & tt<tt.vec[2]){
		scalar=(tt-tt.vec[1])/(tt.vec[2]-tt.vec[1])
	}
	if(tt>tt.vec[3] & tt<tt.vec[4]){
		scalar=(tt.vec[4]-tt)/(tt.vec[4]-tt.vec[3])
	}
	if(tt<=tt.vec[1] | tt>=tt.vec[4]){
		scalar=0
	}
	scalar	
}

#function to calculate temp scalar for growth given temp and temp threasholds;scalar=1.25 for optimal temp, linearly decrease to 0 at upper and lower lethal temps
fun.temp.scalar.growth=function(tt.vec,tt){
	if(tt>=tt.vec[2] & tt<=tt.vec[3]){
		scalar=1.25
	}
	if(tt>tt.vec[1] & tt<tt.vec[2]){
		scalar=1.25*(tt-tt.vec[1])/(tt.vec[2]-tt.vec[1])
	}
	if(tt>tt.vec[3] & tt<tt.vec[4]){
		scalar=1.25*(tt.vec[4]-tt)/(tt.vec[4]-tt.vec[3])
	}
	if(tt<=tt.vec[1] | tt>=tt.vec[4]){
		scalar=0
	}
	scalar	
}

############flow effect
flow.data=read.table("input_flow.txt",header=T, sep='\t',row.names=NULL) #monthly flow time series
flow=flow.data$flow1 #switch between flow scenarios;flow1 is baseline flow setting
surv.flow=2 #m/sec, no effect beblow this flow speed 
max.flow=4 #m/sec, max flow speed assumed in the model, twice the surv.flow
#scalar for nest survival: 1 before surv.flow, linearly decreasing to 0.8 at max.flow
#scalar for fry (<25 mm) mortality: 1 before surv.flow, lnearly increasing to 2 at max.flow

##############chemical effect
chem.on=0 #an indicator to switch on/off chemical effect; 0-off; 1-on

if(chem.on==1){
	chem.data=read.table('input_chemical.txt',header=T,sep='\t',row.names=NULL) #monthly chemcial ng/L time series
	chem=chem.data$chem1 #switch between chem scenarios
	
	chem.level=c(0,3.2,5.3,10.9) #ng/L based on Schwindt et al. 2014
	max.chem=max(chem.level) #maximum chemical level
	fry.surv.chem=c(1,0.38,0.31,0.01) #survival scalar corresponding to chemcial level
	juv.surv.chem=c(1,0.63,0.54,0.01) #survival scalar corresponding to chemcial level
	aduM.surv.chem=c(1,0.26,0.21,0) #survival scalar corresponding to chemcial level
	aduF.surv.chem=c(1,1,0.74,0.68) #survival scalar corresponding to chemcial level
	fecun.chem=c(1,0.65,0.57,0.05) #survival scalar corresponding to chemcial level
	
	fun.chem.scalar=function(y.vec,cc){ #function to calculate chem scalar given scalar at exposure levels
		if(cc>=chem.level[1] & cc<chem.level[2]){
			scalar=y.vec[2]+((y.vec[1]-y.vec[2])*(chem.level[2]-cc)/(chem.level[2]-chem.level[1]))
		}
		if(cc>=chem.level[2] & cc<chem.level[3]){
			scalar=y.vec[3]+((y.vec[2]-y.vec[3])*(chem.level[3]-cc)/(chem.level[3]-chem.level[2]))
		}
		if(cc>=chem.level[3] & cc<=chem.level[4]){
			scalar=y.vec[4]+((y.vec[3]-y.vec[4])*(chem.level[4]-cc)/(chem.level[4]-chem.level[3]))
		}
		scalar
	}
	
	fun.intersex.prob=function(cc){ #function to calculate intersex prob
		prob=0+((1-0)/(1+10^((0.777-cc/0.17)*1.017)))
		prob
	}
	
	fun.deserting.prob=function(cc){ #function to calculate nest deserting prob
		prob=(0+((6-0)/(1+10^((4.492-cc/0.17)*0.911))))*(1/6)
		prob
	}	
}

########################################popn module
########what to track: popn size, PSD, spawner (mature female), recruit (age-0) at monthly scale
stock.size=175 #mm stock length for calculating PSD
quality.size=275 #mm quality size for calculating PSD
quality.low=246 #mm lower bound for quality length
quality.up=279 #mm upper bound for quality length

outcome.popn=matrix(NA,nrow=sim.yr*12,ncol=replicate)
outcome.sp=matrix(NA,nrow=sim.yr*12,ncol=replicate)
outcome.rec=matrix(NA,nrow=sim.yr*12,ncol=replicate)
outcome.psd=matrix(NA,nrow=sim.yr*12,ncol=replicate)
outcome.quality=matrix(NA,nrow=sim.yr*12,ncol=replicate)

#########################################simulation
for(repli in 1:replicate){

##creat table to record individual information, updatae each time step and expand with recruits every year
dd=matrix(NA,nrow=n.initial,ncol=6)
colnames(dd)=c('sex','normalK','mature','spMonth','l','stage') #assing spMonth=0 at begining of each year

dd[,'sex']=rbinom(n.initial, size=1,prob=sex.female) #female-1, male-0
dd[,'normalK']=round(exp(rnorm(n.initial,mean=log(k.mean.normal),sd=logk.sd)),4)
dd[which(dd[,'normalK']<0),'normalK']=0 #no negative growth
dd[,'mature']=0 #assign 0 initially, will be updated at the beginning of each month

##inital size
l.dist=round(n.initial*stru.initial[,'prob']) #num of ind in each size class
l.dist[nrow(stru.initial)]=n.initial-sum(l.dist[1:(nrow(stru.initial)-1)]) #calculate the last size class by -previous size classes,make sure sum up = n.inital after round up
l.resident=NA
for (j in 1:nrow(stru.initial)){
	l.part=round(runif(l.dist[j],min=stru.initial[j,'low'],max=stru.initial[j,'up']),2)
	l.resident=c(l.resident,l.part)
}
l.resident=na.omit(l.resident) #get rid of the NA initially assigned
dd[,'l']=sample(l.resident,size=n.initial,replace=F) #randomize it, not necessary, but better	

########### simulation
for(yy in 1:sim.yr){
	if(nrow(dd)==0){break}
	dd[,'spMonth']=0 #initalize spwning month at the begining of each year, update throughout the year; value for female use as an indicator to track if this female has spawned in this year yet; >0 yes; =0 not yet

for(mm in 1:12){
	if(nrow(dd)==0){break}
	
	#determine maturity:0-immature; 1-mature	
	dd[which(dd[,'l']>mat.size.up & dd[,'mature']==0),'mature']=1
	dd[which(dd[,'l']<mat.size.low & dd[,'mature']==0),'mature']=0
	mature.id=which(dd[,'l']>=mat.size.low & dd[,'l']<=mat.size.up & dd[,'mature']==0)
	if(length(mature.id)>0){
		mat.prob=round((dd[mature.id,'l']-mat.size.low)/(mat.size.up-mat.size.low),3)
		dd[mature.id,'mature']=rbinom(length(mature.id),size=1,prob=mat.prob)
	}

	#assign life stage:1-pre-dispersal fry; 2-post-dispersal fry; 3-juvenile; 4-adults
	dd[which(dd[,'l']<disper.size),'stage']=1 #life stage: 1-pre-dispersal fry; 2-post-dispersal fry; 3-juvenile; 4-adults
	dd[which(dd[,'l']>= disper.size & dd[,'l']<juv.size.low),'stage']=2
	dd[which(dd[,'l']>= juv.size.low & dd[,'mature']==0),'stage']=3
	dd[which(dd[,'mature']==1),'stage']=4
	
	# survival	
	mort.normal=round(m0/dd[,'l'],3)
	
	###temperature effect on mortality
	scalar.temp.mort=rep(1,nrow(dd)) #scalar=1 is no effect
	
	#temp on fry mortality
	scalar.temp.mort[which(dd[,'stage']==1 | dd[,'stage']==2)]=fun.temp.scalar.mort(fry.surv.temp,temp[(yy-1)*12+mm])
		
	#temp on juv mortality
	scalar.temp.mort[which(dd[,'stage']==3)]=fun.temp.scalar.mort(juv.surv.temp,temp[(yy-1)*12+mm])
	
	#temp on adu mortality
	scalar.temp.mort[which(dd[,'stage']==4)]=fun.temp.scalar.mort(adu.surv.temp,temp[(yy-1)*12+mm])

	#density dependence on fry mortality
	scalar.density.mort=rep(1,nrow(dd)) #scalar=1 is no effect
	if(nrow(dd)<(river*max.density*0.5)){
		scalar.density.mort[which(dd[,'stage']==1 | dd[,'stage']==2)]=1
	}
	if(nrow(dd)>(river*max.density)){
		scalar.density.mort[which(dd[,'stage']==1 | dd[,'stage']==2)]=10
	}
	if(nrow(dd)>=(river*max.density*0.5) & nrow(dd)<=(river*max.density)){
		scalar.density.mort[which(dd[,'stage']==1 | dd[,'stage']==2)]=1+((2-1)*(nrow(dd)-river*max.density*0.5)/(river*max.density-river*max.density*0.5))
	}
	
	#flow on fry mortality (l<25 mm)
	scalar.flow.mort=rep(1,nrow(dd)) #scalar=1 is no effect
	if(flow[(yy-1)*12+mm]<surv.flow){
		scalar.flow.mort[which(dd[,'l']<25)]=1
	}
	if(flow[(yy-1)*12+mm]>=surv.flow & flow[(yy-1)*12+mm]<=max.flow){
		scalar.flow.mort[which(dd[,'l']<25)]=1+((2-1)*(flow[(yy-1)*12+mm]-surv.flow)/(max.flow-surv.flow))
	}
	
	#combined mortality before adding chemical effects
	mort=mort.normal*scalar.temp.mort*scalar.density.mort*scalar.flow.mort
	surv.prob=exp((-1)*mort)
	
	#chemical effect on fry, juvenile, adult male and female survival
	if(chem.on==1){ #chemical effect is on
		scalar.chem.surv=rep(1,nrow(dd)) #scalar=1 is no effect
		
		#chemical effect on fry survival
		scalar.chem.surv[which(dd[,'stage']==1 | dd[,'stage']==2)]=fun.chem.scalar(fry.surv.chem,chem[(yy-1)*12+mm])
		
		#chemical effect on juvenile survival
		scalar.chem.surv[which(dd[,'stage']==3)]=fun.chem.scalar(juv.surv.chem,chem[(yy-1)*12+mm])
		
		#chemical effect on adult male survival
		scalar.chem.surv[which(dd[,'stage']==4 & dd[,'sex']==0)]=fun.chem.scalar(aduM.surv.chem,chem[(yy-1)*12+mm])
		
		#chemical effect on adult female survival
		scalar.chem.surv[which(dd[,'stage']==4 & dd[,'sex']==1)]=fun.chem.scalar(aduF.surv.chem,chem[(yy-1)*12+mm])	
	}
	if(chem.on==0){ #chemical effect is off
		scalar.chem.surv=1
	}
	
	#final surv.prob
	surv.prob=surv.prob*scalar.chem.surv
	
	surv=rbinom(nrow(dd),size=1,prob=surv.prob) #1-survive, 0-dead
	surv[which(dd[,'l']>max.size)]=0 #ind > max.size will die
		
	dd=dd[which(surv==1),] #remove dead ind record from table	
	if(length(which(surv==1))==1){dd=t(as.matrix(dd))}	#have to transpose, otherwise R treat it as a vector
	if(nrow(dd)==0){break}
	
	count.survivor=nrow(dd) #record num of remainers for calculate growth later
	outcome.popn[(yy-1)*12+mm,repli]=count.survivor #update popn size after survival
	
	
	#reproduction
	#determine for each male if he is going to spawn in current month based on ddays
	ddays.cal=round(exp(log.ddays.int+log.ddays.slope*log(dd[,'l']/10)),2) #calculate for all ind, but only use those for mature male
	dd[which(dd[,'sex']==0 & dd[,'mature']==1 & dd[,'spMonth']==0 & ddays.cal>ddays.low[(yy-1)*12+mm] & ddays.cal<=ddays.up[(yy-1)*12+mm]),'spMonth']=mm
	
	num.sp.male=length(which(dd[,'sex']==0 & dd[,'spMonth']==mm)) #num of spawning male in current month
	id.sp.female.ava=which(dd[,'sex']==1 & dd[,'mature']==1 & dd[,'spMonth']==0) #available spawning female:mature and not spawn yet (sp.month=0)
	
	if(num.sp.male>0 & length(id.sp.female.ava>0)){
		if(num.sp.male<length(id.sp.female.ava)){
			num.sp.female=num.sp.male
			dd.tem=dd[id.sp.female.ava,]
			dd.tem=cbind(dd.tem,id.sp.female.ava)
			dd.tem=dd.tem[order(dd.tem[,'l'],decreasing=T),]
			id.sp.female=dd.tem[1:num.sp.female,'id.sp.female.ava']
		}else{
			id.sp.female=id.sp.female.ava
		}
		
		dd[id.sp.female,'spMonth']=mm #update female spawning indicator
		outcome.sp[(yy-1)*12+mm,repli]=length(id.sp.female) #record spawner number
						
		fecun=fec.int+fec.slope*dd[id.sp.female,'l'] #num of eggs for each spawning female
		fecun[which(fecun<0)]=0
		
		#chemical effect on fecundity
		if(chem.on==1){
			scalar.chem.fecun=fun.chem.scalar(fecun.chem,chem[(yy-1)*12+mm])
		}
		if(chem.on==0){
			scalar.chem.fecun=1
		}		
			
		#temp effct on egg survival
		scalar.temp.egg=fun.temp.scalar.surv(egg.surv.temp,temp[(yy-1)*12+mm])
					
		#density dependence on egg survival
		if(nrow(dd)<(river*max.density*0.5)){
			scalar.density.egg=1
		}
		if(nrow(dd)>(river*max.density)){
			scalar.density.egg=0
		}
		if(nrow(dd)>=(river*max.density*0.5) & nrow(dd)<=(river*max.density)){
			scalar.density.egg=(river*max.density-nrow(dd))/(river*max.density-river*max.density*0.5)
		}

		egg=fecun*scalar.chem.fecun*fert*egg.surv.normal*scalar.temp.egg*scalar.density.egg #num of fertilized and survived eggs per spawning female (i.e., per nest), a vector
		
		#temp effect on nest survival
		scalar.temp.nest=fun.temp.scalar.surv(nest.surv.temp,temp[(yy-1)*12+mm])
		
		#flow effect on nest survival
		if(flow[(yy-1)*12+mm]<surv.flow){
			scalar.flow.nest=1
		}
		if(flow[(yy-1)*12+mm]>=surv.flow & flow[(yy-1)*12+mm]<=max.flow){
			scalar.flow.nest=0.8+((1-0.8)*(max.flow-flow[(yy-1)*12+mm])/(max.flow-surv.flow))
		}
			
		#chemical effect on nest survival
		if(chem.on==1){ #chemical effect is on
			scalar.chem.nest=1-fun.intersex.prob(chem[(yy-1)*12+mm])*fun.deserting.prob(chem[(yy-1)*12+mm])
		}
		if(chem.on==0){ #chemical effect is off
			scalar.chem.nest=1
		}
			
		nest.surv.prob=nest.surv.normal*scalar.temp.nest*scalar.flow.nest*scalar.chem.nest
		nest.surv=rbinom(length(id.sp.female),size=1,prob=nest.surv.prob) #1-surv, 0-dead
			
		recr=round(sum(egg*nest.surv)) #number of recruits for this month
								
		if(recr>0){
			#density dependence effect
			if((count.survivor+recr)>(river*max.density)){
				recr=river*max.density-count.survivor
			}
				
			dd.new=matrix(NA,nrow=recr,ncol=ncol(dd)) #new recruit table, assign info for new recruits; no need to specify spMonth for new recr,will assign 0 at the beginning of next year
			colnames(dd.new)=colnames(dd)
			
			dd.new[,'sex']=rbinom(recr, size=1,prob=sex.female) #female-1, male-0		
			dd.new[,'normalK']=round(exp(rnorm(recr,mean=log(k.mean.normal),sd=logk.sd)),4)
			dd.new[which(dd.new[,'normalK']<0),'normalK']=0 #no negative growth
			dd.new[,'l']=round(runif(recr,min=swimup.size.low,swimup.size.up),2) #initial length for new recruits
			dd.new[,'mature']=0
			dd.new[,'spMonth']=0
			
			dd=rbind(dd,dd.new) #update dd table with new recruits from this month
				
			outcome.rec[(yy-1)*12+mm,repli]=recr #record recruits number
			outcome.popn[(yy-1)*12+mm,repli]=count.survivor+recr #update popn size with recr
		}
	}
		
	#growth
	##temp on growth rate k
	scalar.temp.growth=fun.temp.scalar.growth(growth.temp,temp[(yy-1)*12+mm])	
	k=dd[1:count.survivor,'normalK']*scalar.temp.growth

	dl=(linf-dd[1:count.survivor,'l'])*(1-exp((-1)*k*dt))
	dd[1:count.survivor,'l']=dd[1:count.survivor,'l']+dl #update length for survivors, not for new recruits
	dd[which(dd[1:count.survivor,'l']>max.size),'l']=max.size
	
	if(length(which(dd[,'l']>=quality.size))>=0 & length(which(dd[,'l']>=stock.size))>0){	
		outcome.psd[(yy-1)*12+mm,repli]=round(length(which(dd[,'l']>=quality.size))/length(which(dd[,'l']>=stock.size))*100,2) #record psd %
	}

	if(length(which(dd[,'l']>=quality.low & dd[,'l']<=quality.up))>=0){	
		outcome.quality[(yy-1)*12+mm,repli]=round(length(which(dd[,'l']>=quality.low & dd[,'l']<=quality.up))/length(dd[,'l']),4) #record quality proportion
	}

} #end mm

write.table(outcome.popn,'popn.txt',quote=F,sep='\t',row.names=F) #report either initial or ending month as outcome this year
write.table(outcome.sp,'sp.txt',quote=F,sep='\t',row.names=F) #need to sum across sp season to get total this year
write.table(outcome.rec,'rec.txt',quote=F,sep='\t',row.names=F) #need to sum across sp season to get total this year
write.table(outcome.psd,'psd.txt',quote=F,sep='\t',row.names=F) #report either initial or ending month as psd this year
write.table(outcome.quality,'quality.txt',quote=F,sep='\t',row.names=F) #report either initial or ending month as outcome this year	
} #end yy
} #end repli

