# pepset

Tryptic peptides similarity among proteomes of multiple organisms

## Get docker image

Either get docker image from docker hub or build docker image locally:

### Pull docker image from docker hub

Check versions here:
```
https://hub.docker.com/repository/docker/mengchen18/pepset
```
Pull the image
```
docker pull mengchen18/pepset:[version number]
```

### Or build docker image locally
```
git clone git@github.com:mengchen18/pepset.git
cd pepset/
docker build . -t pepset
```

## run docker image
```
docker run --rm -v /your/proejct/dir/:/home/shiny/projects -p 4141:3838 pepset
```

## start pepset
open in browser
```
localhost:4141
```
