import os

if __name__ == "__main__":
   f = open(".github/workflows/main.yml","r")
   inp = f.read()
   f.close()

   for line in inp.splitlines() :
       if "replica:" in line :
          adirs = line.replace("replica:","").replace("[","").replace("]","").replace(",","").split()
          rdirs = os.listdir("tests")
          if rdirs != adirs : 
             ValueError("Tests have not been updated. Run the command create_workflow.py in create_workflow directory and commit the changes to .github/actions/main.yml")  
