#!/usr/bin/env bash

# Format all source files in the project.
#
# Set CLANG_FORMAT_CMD environment variable to specify executable used (default: `clang-format`).

set -e

################################################################################
# Configuration

# folders to search
FOLDERS="
    include
    internal
    examples
"

# paths to ignore
EXCLUDE_PATHS="
    internal/ceres/gtest/*
    internal/ceres/gmock/*
    internal/ceres/gmock_gtest_all.cc
    internal/ceres/gmock_main.cc
    internal/ceres/generated/*
    internal/ceres/generated_bundle_adjustment_tests/*
    internal/ceres/schur_eliminator.cc
    internal/ceres/partitioned_matrix_view.cc
    internal/ceres/schur_templates.cc
"

################################################################################
# Implementation

# directory of this script and the repository root
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
REPO_DIR="$SCRIPT_DIR/.."

# set default for CLANG_FORMAT_CMD
CLANG_FORMAT_CMD=${CLANG_FORMAT_CMD:-clang-format}
echo "Formatting with $CLANG_FORMAT_CMD (`$CLANG_FORMAT_CMD --version`)"

# prepare arguments to exclude ignored paths
EXCLUDE_ARGS=""
for p in $EXCLUDE_PATHS; do
    EXCLUDE_ARGS="-not -path */$p $EXCLUDE_ARGS"
done

# for each folder, format header and source dirs
for d in $FOLDERS; do
    d="$REPO_DIR/$d"
    find "$d" \( -name "*.h" -or -name "*.cc" \) $EXCLUDE_ARGS | xargs $CLANG_FORMAT_CMD -verbose -i
done
