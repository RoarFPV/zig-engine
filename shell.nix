
with import <nixpkgs> {};

pkgs.mkShell {
  nativeBuildInputs = with pkgs; [
    zig
    gdb
    SDL2
  ];

  # hardeningDisable = [ "all" ];
}