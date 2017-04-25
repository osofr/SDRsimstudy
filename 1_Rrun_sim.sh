#!/bin/bash
Rfile=$1
# java -Xmx50g -jar ~/h2o-3.8.2.9/h2o.jar -name MyCloud -port 55555 -nthreads 60 > h2o.log &
Rscript --no-save --no-restore --verbose --vanilla --slave $Rfile > $Rfile.out 2>&1
echo "script execution complete"
exit 0



