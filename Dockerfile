FROM mambaorg/micromamba:1.5.10

USER root

RUN apt-get update \
    && apt-get install -y --no-install-recommends procps \
    && rm -rf /var/lib/apt/lists/*

USER $MAMBA_USER

COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml

RUN micromamba install \
    --yes \
    --name base \
    --file /tmp/environment.yml \
    && micromamba clean --all --yes

ENV PATH="/opt/conda/bin:${PATH}"

WORKDIR /workspace

CMD ["bash"]
