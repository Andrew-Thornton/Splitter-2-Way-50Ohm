#!/usr/bin/env bash
set -uo pipefail

EINVAL=22
BASEDIR=$(pwd)/sub/openEMS-Project
VENV_DIR=$(pwd)/venv

DEFAULT_INSTALL_PATH="$HOME/opt/openEMS"

function die {
  printf "%s\n" "A fatal error has occurred!"
  exit 1
}

function build {
  local srcdir; srcdir=$(readlink -f "$1")
  local builddir; builddir=$(readlink -f "$2")
  local njobs="$3"
  local logfile="$4"
  local output="$5"
  local extra_build_arguments=( "${@:6}" )

  cd "$srcdir" || die
  make clean &> /dev/null || true

  if [ -f "$srcdir/bootstrap.sh" ]; then
    echo "bootstrapping $srcdir ... please wait"
    sh ./bootstrap.sh 2>&1 | tee -a "$logfile" >> "$output" || die
  fi

  cd "$builddir" || die

  if [ -f "$srcdir/configure" ]; then
    echo "configuring $srcdir via autotools ... please wait"
    "$srcdir/configure" "${extra_build_arguments[@]}" 2>&1 | tee -a "$logfile" >> "$output" || die
  elif [ -f "$srcdir/CMakeLists.txt" ]; then
    echo "configuring $srcdir via CMake ... please wait"
    cmake "$srcdir" "${extra_build_arguments[@]}" 2>&1 | tee -a "$logfile" >> "$output" || die
  fi

  echo "compiling $srcdir ... please wait"
  make -j"$njobs" 2>&1 | tee -a "$logfile" >> "$output" || die

  cd "$BASEDIR" || die
}

function install {
  local builddir; builddir=$(readlink -f "$1")
  local extra_build_arguments=( "${@:2}" )

  cd "$builddir" || die
  echo "installing $builddir ... please wait"
  make install "${extra_build_arguments[@]}" 2>&1 | tee -a "$LOG_FILE" >> "$STDOUT" || die
  cd "$BASEDIR" || die
}

function setup_venv {
  if [ ! -d "$VENV_DIR" ]; then
    echo "creating python venv at $VENV_DIR"
    python3.12 -m venv "$VENV_DIR" || die
  fi
  source "$VENV_DIR/bin/activate"

  local requirements_file="$(pwd)/requirements.txt"
  if [ -f "$requirements_file" ]; then
    echo "installing python requirements from $requirements_file ... please wait"
    pip install --upgrade pip 2>&1 | tee -a "$LOG_FILE" >> "$STDOUT" || die
    pip install -r "$requirements_file" 2>&1 | tee -a "$LOG_FILE" >> "$STDOUT" || die
  else
    echo "no requirements.txt found at $requirements_file, skipping"
  fi
}

# defaults
STDOUT="/dev/null"
NJOBS=$(python3 -c "import os; print(os.cpu_count())" || nproc || sysctl -n hw.ncpu)
BUILD_HYP2MAT=0
BUILD_CTB=0
BUILD_GUI="YES"
WITH_MPI=0
BUILD_PY_EXT=1
BUILD_TINYXML=0
INSTALL_PATH="$DEFAULT_INSTALL_PATH"
PYTHON_ARGS=()
LOG_FILE="$BASEDIR/build_$(date +%Y%m%d_%H%M%S).log"

if [[ "$OSTYPE" == "darwin"* ]]; then
  BUILD_TINYXML=1
fi

echo "setting install path to: $INSTALL_PATH"
echo "logging build output to: $LOG_FILE"

setup_venv

TMPDIR=$(mktemp -d)
mkdir -p "$INSTALL_PATH"

# build TinyXML
if [ "$BUILD_TINYXML" -eq 1 ]; then
  echo "downloading and building custom TinyXML in tmp dir: $TMPDIR"
  mkdir -p ./downloads
  ./scripts/build_tinyxml.sh --build-dir "$TMPDIR" --install-dir "$INSTALL_PATH" \
    2>&1 | tee -a "$LOG_FILE" >> "$STDOUT" || die
fi

# build openEMS Project
build "$BASEDIR" "$TMPDIR" "$NJOBS" "$LOG_FILE" "$STDOUT" \
      "-DBUILD_APPCSXCAD=$BUILD_GUI" \
      "-DCMAKE_INSTALL_PREFIX=$INSTALL_PATH" \
      "-DWITH_MPI=$WITH_MPI"

# hyp2mat
if [ $BUILD_HYP2MAT -eq 1 ]; then
  mkdir -p "$TMPDIR/hyp2mat"
  build hyp2mat "$TMPDIR/hyp2mat" "$NJOBS" "$LOG_FILE" "$STDOUT" "--prefix=$INSTALL_PATH"
  install "$TMPDIR/hyp2mat"
fi

# circuit toolbox
if [ $BUILD_CTB -eq 1 ]; then
  install CTB "PREFIX=$INSTALL_PATH"
fi

# python extension build
echo "Building python modules ... please wait"
./scripts/build_python.sh \
  --cpp-install-dir "$INSTALL_PATH" \
  ${PYTHON_ARGS[@]+"${PYTHON_ARGS[@]}"} \
  2>&1 | tee -a "$LOG_FILE" >> "$STDOUT" || die

cd "$BASEDIR" || die

echo "build successful, cleaning up tmp dir ..."
rm -rf "$TMPDIR"

echo ""
echo "% the python venv used for this build can be activated any time with:"
echo "source ./venv/bin/activate"
echo ""
