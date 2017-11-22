# Docker image workflow
- First run `make <build_flavour>_<target>_docker` on the local environment ICC toolchain, since it can't be redistributed inside a 'build' Docker image
- Then use the CLI to build the Docker image of the sim, using the fresh binaries it depends on, and either specifying a tag or keeping the ID of the generated image (which makes make-based rebuilds awkward, for now)
- Use the previous tag or ID to start a container, inside which `/app/bin/infoli.x.docker` can be invoked directly, through `sh`, or later on a work driver process.
