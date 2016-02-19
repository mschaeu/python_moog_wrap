#     makefile for MOOG with all of the common block assignments;

#     here are the object files
EXEC   := moog.so
SRC    := Abfind.f Abpop.f Batom.f Begin.f Binary.f \
	Blankstring.f Blends.f Bmolec.f \
	Calmod.f Cdcalc.f Chabund.f Cog.f Cogsyn.f \
	Curve.f Damping.f Discov.f \
	Doflux.f Eqlib.f Estim.f Ewfind.f \
	Ewweighted.f Fakeline.f Finish.f \
	Gammabark.f Getasci.f Getcount.f Getnum.f Getsyns.f \
	Gridplo.f Gridsyn.f Infile.f Inlines.f Inmodel.f \
	Inmodel_python.f Invert.f \
	Jexpint.f Lineinfo.f Lineabund.f Linlimit.f \
	Minimax.f Molquery.f Moog_lib.f Mydriver.f \
	Nansi.f Nearly.f Number.f Obshead.f\
	Oneline.f Opaccouls.f OpacHelium.f OpacHydrogen.f \
	Opacit.f Opacmetals.f Opacscat.f Params.f Partfn.f \
	Partnew.f Plotit.f \
	Prinfo.f Putasci.f Readobs.f \
	Rinteg.f Smooth.f Stats.f Sunder.f Synpop.f Synspec.f \
	Synth.f Tablepop.f Taukap.f Total.f Trudamp.f Ucalc.f \
	Vargauss.f Vmacro.f Voigt.f Wavecalc.f Weedout.f Writenumber.f	

#     here are the common files
COMMON =  Atmos.com Dummy.com Equivs.com Factor.com Kappa.com Linex.com \
	Mol.com Multistar.com Obspars.com Plotval.com Pstuff.com \
	Quants.com Multimod.com Dampdat.com
	

FC = gfortran -w

$(EXEC): $(SRC)
	$(FC) -ffree-form -shared $^ -o $@


$(SRC): $(COMMON)

# Phony targets
.PHONY: clobber clean neat echo run
clean: neat
	$(RM) $(EXEC) 
neat:
	$(RM) *~ .*~
echo:
	@echo $(OBJ)
