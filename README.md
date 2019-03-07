# fastslam_web
Experimental fastslam algorithm webapp.

## Installation

To deploy locally:

1. Install [Docker](https://www.docker.com/get-started) on your machine.

2. Clone the repository

```
git clone https://github.com/pptacher/fastslam_web.git
```

3. Go to the root of the repository

```
cd fastslam_web
```

4. Build the docker image

```
docker build --tag=*fastslam* .
```

5. Run it

```
docker -p 3000:3000 run *fastslam*
```

6. Navigate to [http://localhost:3000](http://localhost:3000). If you use Docker through Linux virtual machine (MacOS), you may have to access the web app with something like [192.168.99.100:3000](192.168.99.100:3000).
