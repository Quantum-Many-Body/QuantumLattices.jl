using CompileBot

botconfig = BotConfig(
  "QuantumLattices";                   # package name (the one this configuration lives in)
  yml_path = "SnoopCompile.yml"        # parse `os` and `version` from `SnoopCompile.yml`
)

snoop_bot(
  botconfig,
  "$(@__DIR__)/example_script.jl",
)