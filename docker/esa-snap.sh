#!/bin/sh

wget -q -O /src/snap_install.sh "http://step.esa.int/downloads/9.0/installers/esa-snap_sentinel_unix_9_0_0.sh"
sh /src/snap_install.sh -q -varfile /src/docker/esa-snap.varfile

# https://senbox.atlassian.net/wiki/spaces/SNAP/pages/30539785/Update+SNAP+from+the+command+line
/usr/local/snap/bin/snap --nosplash --nogui --modules --update-all 2>&1 | while read -r line; do
    echo "$line"
    [ "$line" = "updates=0" ] && sleep 2 && pkill -TERM -f "snap/jre/bin/java"
done

rm -rf /src/snap_install.sh
