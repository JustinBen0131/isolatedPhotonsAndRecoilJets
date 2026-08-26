# Workflows

This directory will define the inspectable stage graph for production,
training, validation, offline analysis, and plotting.

The graph coordinates directly callable commands; it does not require a
workflow server. Local examples and batch execution consume the same stage
contracts and differ only in resource and input-location configuration.
