import os

if __name__ == "__main__":
   f = open("base_workflow.yml","r")
   inp = f.read()
   f.close()

   of = open("../.github/workflows/main.yml","w+")
   for line in inp.splitlines() :
       if "replica:" in line :
          alltests = os.listdir("../tests")
          of.write( line + " [" + alltests[0] )
          for i in range(1,len(alltests)) : of.write(", " + alltests[i])
          of.write("]\n")
       else : of.write( line + "\n")
   of.close()
