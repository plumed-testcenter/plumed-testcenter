import subprocess

class mdcode :
   def getName( self ) :
       return "quantum_expresso"

   def setParams( self ) :
       params = {
         "temperature": 1.0,
         "tstep": 20,
         "friction": 1.0
       }
       return params

   def runMD( self, mdparams ) : 
       of = open("md.in","w+")
       of.write(" &control \n")
       of.write("    calculation='md' \n")
       of.write("    pseudo_dir='./' \n")
       of.write("    dt=" + str(mdparams["tstep"]) + ", \n")
       of.write("    nstep=" + str(mdparams["nsteps"]) + " \n")
       of.write(" / \n")
       of.write(" &system \n")
       of.write("    ibrav= 2, celldm(1)=10.18, nat=  2, ntyp= 1, \n")
       of.write("    ecutwfc = 8.0, nosym=.true. \n")
       of.write(" / \n")
       of.write(" &electrons \n")
       of.write("    conv_thr =  1.0e-8 \n")
       of.write("    mixing_beta = 0.7 \n")
       of.write(" / \n")
       of.write(" &ions \n")
       of.write(" / \n")
       of.write("ATOMIC_SPECIES \n")
       of.write(" Si  28.086  Si.pz-vbc.UPF \n")
       of.write("ATOMIC_POSITIONS {alat} \n")
       of.write(" Si -0.123 -0.123 -0.123 \n")
       of.write(" Si  0.123  0.123  0.123 \n")
       of.write("K_POINTS {automatic} \n")
       of.write(" 1 1 1 0 0 0 \n")
       of.close()
