#!/bin/bash
java -$1 -jar  $(dirname "$0")/MSGFPLUS_v20220418/MSGFPlus.jar ${@:2}
