GF = gfortran
object = addopt.o phases-n2lo1opt2.o N2LOFUN1opt.o n2lo500add.o SAphase.o
module = addterms.mod const.mod kqxy.mod
# generate two execute files
# SA.exe
SA.exe: $(object) 
	$(GF) -o SA.exe $(object)
SAphase.o: SAphase.f phase.mod
	$(GF) -c SAphase.f
n2lo500add.o: n2lo500add.f
	$(GF) -c n2lo500add.f
N2LOFUN1opt.o: N2LOFUN1opt.f
	$(GF) -c N2LOFUN1opt.f
phases-n2lo1opt2.o phase.mod: phases-n2lo1opt2.f
	$(GF) -c phases-n2lo1opt2.f
addopt.o $(module): addopt.f
	$(GF) -c addopt.f
.PHONY :clean
clean:
	del *.exe *.o *.mod