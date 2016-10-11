import os

print "pythonpath in script"
print os.environ["PYTHONPATH"]
print "Importing logD"
import predictlogD
print dir(predictlogD)
print "Importing Peff"
import predictPeff
print dir(predictPeff)
