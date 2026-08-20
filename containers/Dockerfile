# One immutable image containing the pipeline and the four external tools, for a machine
# that has no network at all.
#
# Build once on a machine that does have network, carry the result across, run it there
# forever. The point is not convenience: it is that the image digest pins the tools, their
# shared libraries and the operating-system layer together, so "which version produced this
# report" has one answer that survives the machine being rebuilt.
#
#     docker build -f containers/Dockerfile -t ontseq:0.3.0 .
#     docker save ontseq:0.3.0 | gzip > ontseq-0.3.0.tar.gz
#     sha256sum ontseq-0.3.0.tar.gz          # carry this number separately from the medium
#
# On the isolated machine:
#
#     sha256sum -c ontseq-0.3.0.tar.gz.sha256
#     gunzip -c ontseq-0.3.0.tar.gz | docker load
#
# Apptainer, which most clinical compute environments prefer, reads the same artifact:
#
#     apptainer build ontseq-0.3.0.sif docker-archive://ontseq-0.3.0.tar
#
# There is deliberately no second Apptainer definition file. Two build recipes for one image
# are two things that drift apart, and the first sign of drift is a report whose provenance
# names tool versions that were never in the container that produced it.

FROM mambaorg/micromamba:2.0.5

# The same environment file the CI smoke job installs from. Single source of truth for the
# tool versions, so the container and CI cannot end up pinned differently.
COPY --chown=$MAMBA_USER:$MAMBA_USER workflow/envs/aligned_bam.yaml /tmp/aligned_bam.yaml
RUN micromamba install -y -n base -f /tmp/aligned_bam.yaml \
    && micromamba clean --all --yes

ARG MAMBA_DOCKERFILE_ACTIVATE=1

COPY --chown=$MAMBA_USER:$MAMBA_USER . /opt/ontseq-platform

# --no-deps is deliberate. Every runtime dependency is already pinned in the environment
# file and resolved by conda. Letting pip resolve them a second time would install a
# different set of versions on top, and the provenance record would then name versions that
# are not the ones being executed.
RUN pip install --no-deps /opt/ontseq-platform

# So that the default relative policy paths (configs/qc/defaults.yaml and friends) resolve.
# Mount your own data somewhere absolute and refer to it by absolute path.
WORKDIR /opt/ontseq-platform

ENTRYPOINT ["/usr/local/bin/_entrypoint.sh", "ontseq"]
CMD ["--help"]
