How to install: https://hub.docker.com/r/ugrbioinfo/srnatoolbox

Standalone version
```bash
git clone https://github.com/bioinfoUGR/sRNAtoolbox.git
```

To use:
```bash
# sudo docker stop sRNAtoolbox (if wish to stop)
sudo docker start sRNAtoolbox 
sudo docker exec -term=SCREEN -it sRNAtoolbox /bin/bash

# sudo docker run --hostname sRNAtoolbox --name sRNAtoolbox --user srna --workdir $PWD -it ugrbioinfo/srnatoolbox:latest /bin/bash
```