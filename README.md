# pepset

Tryptic peptides similarity among proteomes of multiple organisms

## build docker
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
