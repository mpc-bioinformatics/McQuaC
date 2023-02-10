#!/bin/bash
java -$1 -cp $(dirname "$0")/MSGFPLUS_v20220418/MSGFPlus.jar edu.ucsd.msjava.msdbsearch.BuildSA ${@:2}
