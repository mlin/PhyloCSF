OPAM_DEPENDS="batteries gsl ocaml+twt ounit should"
 
sudo add-apt-repository -y ppa:avsm
sudo apt-get update -qq
sudo apt-get install -qq libgsl0-dev ocaml ocaml-native-compilers camlp4-extra opam
export OPAMYES=1
echo OCaml version
ocaml -version
echo OPAM versions
opam --version

opam init 
opam install ${OPAM_DEPENDS}
eval `opam config env`
bash -c "while :; do date; sleep 60; done" &
killme=$!
make test
kill $killme || true
