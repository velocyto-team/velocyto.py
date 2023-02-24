from python:3.10-slim

ENV DEBIAN_FRONTEND=noninteractive
RUN apt update --fix-missing && \
   apt install --without-recommends -y -f --no-gui \
      build-essential \
      curl \
      git-core \
      libglu1-mesa-dev \
      pandoc \
      pandoc-citeproc \
      libcgal-dev \
      libgit2-dev \
      libgsl-dev \
      libfreetype6-dev \
      libopenblas-dev \
      libssh2-1-dev \
      libxkbcommon-x11-0 \
      rclone \
      ssh \
      vim && \
   cp /usr/share/zoneinfo/America/Chicago /etc/localtime && \
   aptitude clean && \
   rm -rf /tmp/downloaded_packages/* && \
   rm -rf /var/lib/apt/lists/*