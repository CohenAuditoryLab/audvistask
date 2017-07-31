#!/bin/bash
echo "Backing up data ..."
rsync -auv /Users/briannakarpowicz/Documents/Cohen\ Lab/Auditory-Visual\ Task/ jaejin@130.91.169.205:/Volumes/KittyPride/VisualTones/
echo "Done!"
exit
