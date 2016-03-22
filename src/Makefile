CC = g++
CXXFLAGS = -DHAVE_CONFIG_H -O2 -fomit-frame-pointer -funroll-loops
ARCH_AVX = -O2 -mavx -D_TOM_AVX
ARCH_SSE = -O2 -msse3 -D_TOM_SSE3
PAR_OMP  = -fopenmp
ASM_DBG  = -D_ASM_DEBUG
OBJS 	 = dppdiv.o Alignment.o MbEigensystem.o MbMath.o MbRandom.o MbTransitionMatrix.o Mcmc.o Parameter.o Parameter_basefreq.o Parameter_exchangeability.o Parameter_fossilgraph.o Parameter_fossilrangegraph.o Parameter_origin.o Parameter_rate.o Parameter_shape.o Parameter_tree.o Parameter_cphyperp.o Parameter_treescale.o Parameter_speciaton.o Parameter_expcalib.o Calibration.o Model.o
RM 	 = rm -f
PROF	 = -pg
DEBUG    = -DDEBUG -g -O2 -fomit-frame-pointer -funroll-loops


all: dppdiv-seq dppdiv-seq-avx dppdiv-seq-sse dppdiv-par dppdiv-par-sse dppdiv-par-avx
asm: asm-seq asm-seq-avx asm-seq-sse
prof: dppdiv-prof-seq
debug: dppdiv-debug

dppdiv-seq: $(OBJS) Model_likelihood-seq.o
	$(CC) -o $@ $+ -g

dppdiv-seq-avx: $(OBJS) Model_likelihood-seq-avx.o
	$(CC) -O2 -o $@ $(ARCH_AVX) $+

dppdiv-seq-sse: $(OBJS) Model_likelihood-seq-sse.o
	$(CC) -o $@ $(ARCH_SSE) $+

dppdiv-par: $(OBJS) Model_likelihood-par.o
	$(CC) -o $@ $(PAR_OMP) $+

dppdiv-par-sse: $(OBJS) Model_likelihood-par-sse.o
	$(CC) -o $@ $(PAR_OMP) $(ARCH_SSE) $+

dppdiv-par-avx: $(OBJS) Model_likelihood-par-avx.o
	$(CC) -o $@ $(PAR_OMP) $(ARCH_AVX) $+

asm-seq: Model_likelihood.cpp
	$(CC) -S -O2 -o dppdiv-seq.s $(ASM_DBG) $+

asm-seq-avx: Model_likelihood.cpp
	$(CC) -S -o dppdiv-seq-avx.s $(ARCH_AVX) $(ASM_DBG) $+

asm-seq-sse: Model_likelihood.cpp
	$(CC) -S -o dppdiv-seq-sse.s $(ARCH_SSE) $(ASM_DBG) $+

dppdiv-prof-seq: $(OBJS) Model_likelihood-seq.o
	$(CC) $(PROF) -o $@ $+

dppdiv-debug: Alignment.cpp Calibration.cpp dppdiv.cpp MbEigensystem.cpp MbMath.cpp MbRandom.cpp MbTransitionMatrix.cpp Mcmc.cpp Model-old.cpp Parameter_basefreq.cpp Parameter_cphyperp.cpp Parameter.cpp Parameter_exchangeability.cpp Parameter_fossilgraph.cpp Parameter_fossilrangegraph.cpp Parameter_origin.cpp Parameter_expcalib.cpp Parameter_rate.cpp Parameter_shape.cpp Parameter_speciaton.cpp Parameter_tree.cpp Parameter_treescale.cpp
	$(CC) $(DEBUG) -o $@ $+


dppdiv.o: dppdiv.cpp
Alignment.o: Alignment.cpp
MbEigensystem.o: MbEigensystem.cpp
MbMath.o: MbMath.cpp
MbRandom.o: MbRandom.cpp
MbTransitionMatrix.o: MbTransitionMatrix.cpp
Mcmc.o: Mcmc.cpp
Model.o: Model.cpp
Model_likelihood-seq.o: Model_likelihood.cpp
	$(CC) -c -o $@ $(CXXFLAGS) $+
Model_likelihood-par.o: Model_likelihood.cpp
	$(CC) -c -o $@ $(CXXFLAGS) $(PAR_OMP) $+
Model_likelihood-par-sse.o: Model_likelihood.cpp
	$(CC) -c -o $@ $(CXXFLAGS) $(PAR_OMP) $(ARCH_SSE) $+
Model_likelihood-par-avx.o: Model_likelihood.cpp
	$(CC) -c -o $@ $(CXXFLAGS) $(PAR_OMP) $(ARCH_AVX) $+
Model_likelihood-seq-sse.o: Model_likelihood.cpp
	$(CC) -c -o $@ $(CXXFLAGS) $(ARCH_SSE) $+
Model_likelihood-seq-avx.o: Model_likelihood.cpp
	$(CC) -c -o $@ $(CXXFLAGS) $(ARCH_AVX) $+
Parameter.o: Parameter.cpp
Parameter_basefreq.o: Parameter_basefreq.cpp
Parameter_exchangeability.o: Parameter_exchangeability.cpp
Parameter_fossilgraph.o: Parameter_fossilgraph.cpp
Parameter_fossilrangegraph.o: Parameter_fossilrangegraph.cpp
Parameter_origin.o: Parameter_origin.cpp
Parameter_rate.o: Parameter_rate.cpp
Parameter_shape.o: Parameter_shape.cpp
Parameter_tree.o: Parameter_tree.cpp
Parameter_cphyperp.o: Parameter_cphyperp.cpp
Parameter_speciaton.o: Parameter_speciaton.cpp
Parameter_treescale.o: Parameter_treescale.cpp
Parameter_expcalib.o: Parameter_expcalib.cpp
Calibration.o: Calibration.cpp

clean:
	$(RM) *.o dppdiv-seq dppdiv-seq-avx dppdiv-seq-sse dppdiv-par dppdiv-par-sse dppdiv-par-avx dppdiv-seq.s dppdiv-seq-avx.s dppdiv-seq-sse.s
