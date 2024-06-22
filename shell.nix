# -DZIG_STATIC

with import <nixpkgs> { };

pkgs.mkShell {
  nativeBuildInputs = with pkgs; [
    pkg-config
    autoPatchelfHook
    installShellFiles
    #rives
    # - machine-emulator
    boost
    protobuf

    wget

    cmake
    gdb
    libxml2
    ninja
    qemu
    wasmtime
    zlib
    zstd
    libgcc.lib
    SDL2
  ] ++ (with llvmPackages_18; [
    clang
    clang-unwrapped
    lld
    llvm
    libclang
  ]);

  hardeningDisable = [ "all" ];


  shellHook = ''
    export PATH=~/zig/0.13.0/files:$PATH
  '';

}
