# fastslam cpp implementation

Implementation of the fastslam algorithm. You can execute natively after installing locally all the dependencies, or execute it in docker container.

## Building the Docker image and executing in container

1. Install [Docker](https://www.docker.com/get-started) on your machine.

2. Clone the repository

```
git clone https://github.com/pptacher/probabilistic_robotics.git
```

3. Go to the source folder

```
cd probabilistic_robotics/ch13_the_fastslam_algorithm/src/cpp
```

4. Build the docker image

```
docker build --tag=fastslam .
```

5. Run the container; you can specify a number of particles (default to 100).

```
docker run fastslam 150
```

You can access the container id with
```
$ docker container ls
CONTAINER ID        IMAGE               COMMAND             CREATED              STATUS              PORTS                    NAMES
ca6183eab770        fastslam            "./fslam.bin"       About a minute ago   Up About a minute                            unruffled_hellman
```
copy back the output locally to your machine
```
$ docker cp e31da5e0537b:/app/data/poses.dat .
```

pause the container
```
docker pause ca6183eab770
```

stop the container
```
docker stop ca6183eab770
```

following command will remove all stopped containers, all dangling images, and all unused networks:
```
docker system prune
```

See also the [Docker cheat sheet](https://www.docker.com/sites/default/files/Docker_CheatSheet_08.09.2016_0.pdf).
