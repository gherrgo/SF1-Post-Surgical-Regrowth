mixed PostoperativeVolume i.DNAmethylationclass i.Yearsafter1postoperativesca i.DNAmethylationclass##i.Yearsafter1postoperativesca PreoperativeTumorVolume InvasiveGrowth Reintervention || Recordid:
predict res,res
qnorm res 
margins i.DNAmethylationclass, at(Yearsafter1postoperativesca=(1(1)10))
marginsplot
drop if DNAmethylationclass == 1

by Recordid, sort : egen float Maxmonths = max(Monthsafter1postoperativesc)
by Recordid, sort : egen float Maxyears = max(Yearsafter1postoperativesca)
by Recordid, sort : egen float Reinterventionagain = max(Reintervention)
duplicates drop Recordid, force

gen byte InvasiveGrowth_byte = .
replace InvasiveGrowth_byte = 1 if InvasiveGrowth == "1"
replace InvasiveGrowth_byte = 0 if InvasiveGrowth == "0"

summarize Maxmonths, detail
mean PreoperativeTumorVolume, over(DNAmethylationclass)


mixed PostoperativeVolume i.Kcluster_5 i.Yearsafter1postoperativesca i.Kcluster_5##i.Yearsafter1postoperativesca PreoperativeTumorVolume InvasiveGrowth Reintervention || Recordid:
predict res,res
qnorm res 
margins i.Kcluster_5, at(Yearsafter1postoperativesca=(1(1)10))
marginsplot
drop if Kcluster_5 == 4

generate Nogrowth = (PostoperativeVolume > PostoperativeTumorVolume)
stset Monthsafter1postoperativesc, id(Recordid) failure(Nogrowth==1)
stsum, by(Kcluster_5)
sts graph, by(Kcluster_5)
sts graph, by(Kcluster_5) risktable risktable(24 48 72 96 120 144 168, failevents lastfailure) xlabel(0(24)156)

stcox i.Kcluster_5

drop if Kcluster_5 == 4
