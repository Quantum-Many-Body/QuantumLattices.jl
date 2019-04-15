using Test
using QuantumLattices.Essentials.Extensions
using QuantumLattices.Essentials.FockPackage: σ⁰,σˣ,σʸ,σᶻ,σ⁺,σ⁻
using QuantumLattices.Essentials.SpinPackage: Heisenberg,Ising,Gamma,Sˣ,Sʸ,Sᶻ

@testset "fockcoupling" begin
    fc=fc"1.0 sl(1,1,1,1)⊗ob(α,α,β,β)⊗sp(σ,γ,σ,γ)⊗ph(2,1,2,1)@(1,1,2,2) with (α<β,σ≠γ)"
    ob,sp=fc.obsubscripts[1].identifier,fc.spsubscripts[1].identifier
    @test repr(fc)=="1.0 sl(1,1,1,1)⊗ob(α,α,β,β)⊗sp(σ,γ,σ,γ)⊗ph(2,1,2,1)@(1,1,2,2) with ($ob,$sp) && ($ob,$sp) && ($ob,$sp) && ($ob,$sp)"

    fc=fc"1.0 sl(1,1,1,1)⊗ob(α,α,β,β)⊗sp(σ,γ,σ,γ)⊗ph(2,1,2,1)@(1,1,2,2)"
    ob,sp=fc.obsubscripts[1].identifier,fc.spsubscripts[1].identifier
    @test repr(fc)=="1.0 sl(1,1,1,1)⊗ob(α,α,β,β)⊗sp(σ,γ,σ,γ)⊗ph(2,1,2,1)@(1,1,2,2) with ($ob,$sp) && ($ob,$sp) && ($ob,$sp) && ($ob,$sp)"

    fc=fc"1.0 sl(1,1,1,1)⊗ob(α,α,β,β)⊗sp(σ,γ,σ,γ)⊗ph(2,1,2,1)"
    ob,sp=fc.obsubscripts[1].identifier,fc.spsubscripts[1].identifier
    @test repr(fc)=="1.0 sl(1,1,1,1)⊗ob(α,α,β,β)⊗sp(σ,γ,σ,γ)⊗ph(2,1,2,1) with ($ob,$sp) && ($ob,$sp) && ($ob,$sp) && ($ob,$sp)"

    fc=fc"1.0 sl(1,1,1,1)⊗ob(α,α,β,β)⊗sp(σ,γ,σ,γ)⊗ph(2,1,2,1) with (α<β,σ≠γ)"
    ob,sp=fc.obsubscripts[1].identifier,fc.spsubscripts[1].identifier
    @test repr(fc)=="1.0 sl(1,1,1,1)⊗ob(α,α,β,β)⊗sp(σ,γ,σ,γ)⊗ph(2,1,2,1) with ($ob,$sp) && ($ob,$sp) && ($ob,$sp) && ($ob,$sp)"

    fc=fc"1.0 sl(1,1,1,1)⊗ob(α,α,β,β)⊗sp(σ,γ,σ,γ)⊗ph(2,1,2,1) with α<β,σ≠γ"
    ob,sp=fc.obsubscripts[1].identifier,fc.spsubscripts[1].identifier
    @test repr(fc)=="1.0 sl(1,1,1,1)⊗ob(α,α,β,β)⊗sp(σ,γ,σ,γ)⊗ph(2,1,2,1) with ($ob,$sp) && ($ob,$sp) && ($ob,$sp) && ($ob,$sp)"

    fc=fc"1.0 sl(1,1,1,1)⊗ob(α,α,β,β)⊗sp(σ,γ,σ,γ)⊗ph(2,1,2,1) with (α<β,*)"
    ob,sp=fc.obsubscripts[1].identifier,fc.spsubscripts[1].identifier
    @test repr(fc)=="1.0 sl(1,1,1,1)⊗ob(α,α,β,β)⊗sp(σ,γ,σ,γ)⊗ph(2,1,2,1) with ($ob,$sp) && ($ob,$sp) && ($ob,$sp) && ($ob,$sp)"

    fc=fc"1.0 sl(1,1,1,1)⊗ob(α,α,β,β)⊗ph(2,1,2,1) with α<β"
    ob,sp=fc.obsubscripts[1].identifier,fc.spsubscripts[1].identifier
    @test repr(fc)=="1.0 sl(1,1,1,1)⊗ob(α,α,β,β)⊗ph(2,1,2,1) with ($ob,$sp) && ($ob,$sp) && ($ob,$sp) && ($ob,$sp)"

    fc=fc"1.0 sl(1,1,1,1)⊗ob(α,α,β,β)⊗sp(2,1,2,1) with α<β"
    ob,sp=fc.obsubscripts[1].identifier,fc.spsubscripts[1].identifier
    @test repr(fc)=="1.0 sl(1,1,1,1)⊗ob(α,α,β,β)⊗sp(2,1,2,1) with ($ob,$sp) && ($ob,$sp) && ($ob,$sp) && ($ob,$sp)"

    fc=fc"1.0 sl(1,1,1,1)⊗ob(1,1,1,1)⊗ph(2,1,2,1)"
    @test repr(fc)=="1.0 sl(1,1,1,1)⊗ob(1,1,1,1)⊗ph(2,1,2,1)"

    fc=fc"1.0 sl(1,1,1,1)⊗ph(2,1,2,1)"
    @test repr(fc)=="1.0 sl(1,1,1,1)⊗ph(2,1,2,1)"

    fc=fc"1.0 @(1,1,1,1)"
    @test repr(fc)=="1.0 @(1,1,1,1)"

    fc=fc"1.0im {2}"
    @test repr(fc)=="1.0im {2}"
end

@testset "fockcouplings" begin
    @test σ⁰"sp"==σ⁰("sp") && σ⁰"sp@(1,2)"==σ⁰("sp",centers=(1,2))
    @test σ⁰"ob"==σ⁰("ob") && σ⁰"ob@(1,2)"==σ⁰("ob",centers=(1,2))
    @test σ⁰"sl"==σ⁰("sl") && σ⁰"sl@(1,2)"==σ⁰("sl",centers=(1,2))
    @test σ⁰"ph"==σ⁰("ph") && σ⁰"ph@(1,2)"==σ⁰("ph",centers=(1,2))

    @test σˣ"sp"==σˣ("sp") && σˣ"sp@(1,2)"==σˣ("sp",centers=(1,2))
    @test σˣ"ob"==σˣ("ob") && σˣ"ob@(1,2)"==σˣ("ob",centers=(1,2))
    @test σˣ"sl"==σˣ("sl") && σˣ"sl@(1,2)"==σˣ("sl",centers=(1,2))
    @test σˣ"ph"==σˣ("ph") && σˣ"ph@(1,2)"==σˣ("ph",centers=(1,2))

    @test σʸ"sp"==σʸ("sp") && σʸ"sp@(1,2)"==σʸ("sp",centers=(1,2))
    @test σʸ"ob"==σʸ("ob") && σʸ"ob@(1,2)"==σʸ("ob",centers=(1,2))
    @test σʸ"sl"==σʸ("sl") && σʸ"sl@(1,2)"==σʸ("sl",centers=(1,2))
    @test σʸ"ph"==σʸ("ph") && σʸ"ph@(1,2)"==σʸ("ph",centers=(1,2))

    @test σᶻ"sp"==σᶻ("sp") && σᶻ"sp@(1,2)"==σᶻ("sp",centers=(1,2))
    @test σᶻ"ob"==σᶻ("ob") && σᶻ"ob@(1,2)"==σᶻ("ob",centers=(1,2))
    @test σᶻ"sl"==σᶻ("sl") && σᶻ"sl@(1,2)"==σᶻ("sl",centers=(1,2))
    @test σᶻ"ph"==σᶻ("ph") && σᶻ"ph@(1,2)"==σᶻ("ph",centers=(1,2))

    @test σ⁺"sp"==σ⁺("sp") && σ⁺"sp@(1,2)"==σ⁺("sp",centers=(1,2))
    @test σ⁺"ob"==σ⁺("ob") && σ⁺"ob@(1,2)"==σ⁺("ob",centers=(1,2))
    @test σ⁺"sl"==σ⁺("sl") && σ⁺"sl@(1,2)"==σ⁺("sl",centers=(1,2))
    @test σ⁺"ph"==σ⁺("ph") && σ⁺"ph@(1,2)"==σ⁺("ph",centers=(1,2))

    @test σ⁻"sp"==σ⁻("sp") && σ⁻"sp@(1,2)"==σ⁻("sp",centers=(1,2))
    @test σ⁻"ob"==σ⁻("ob") && σ⁻"ob@(1,2)"==σ⁻("ob",centers=(1,2))
    @test σ⁻"sl"==σ⁻("sl") && σ⁻"sl@(1,2)"==σ⁻("sl",centers=(1,2))
    @test σ⁻"ph"==σ⁻("ph") && σ⁻"ph@(1,2)"==σ⁻("ph",centers=(1,2))
end

@testset "spincoupling" begin
    sc=sc"1.0 S+S- sl(1,1)⊗ob(α,β)@(1,2) with α<β"
    ob=sc.subscripts[1].identifier
    @test repr(sc)=="1.0 S+S- sl(1,1)⊗ob(α,β)@(1,2) with $ob && $ob"

    sc=sc"1.0 S+S- sl(1,1)⊗ob(α,β)@(1,2)"
    ob=sc.subscripts[1].identifier
    @test repr(sc)=="1.0 S+S- sl(1,1)⊗ob(α,β)@(1,2) with $ob && $ob"

    sc=sc"1.0 S+S- sl(1,1)⊗ob(α,β)"
    ob=sc.subscripts[1].identifier
    @test repr(sc)=="1.0 S+S- sl(1,1)⊗ob(α,β) with $ob && $ob"

    sc=sc"1.0 S+S- sl(1,1)⊗ob(1,2)@(1,2)"
    @test repr(sc)=="1.0 S+S- sl(1,1)⊗ob(1,2)@(1,2)"

    sc=sc"1.0 S+S- sl(1,1)@(1,2)"
    @test repr(sc)=="1.0 S+S- sl(1,1)@(1,2)"

    sc=sc"1.0 S+S- @(1,2)"
    @test repr(sc)=="1.0 S+S- @(1,2)"

    sc=sc"1.0 S+S-"
    @test repr(sc)=="1.0 S+S-"
end

@testset "spincouplings" begin
    @test heisenberg"sl(1,1)⊗ob(1,3)@(1,2)"==Heisenberg(centers=(1,2),atoms=(1,1),orbitals=(1,3))
    @test heisenberg"sl(1,1)⊗ob(1,3)"==Heisenberg(atoms=(1,1),orbitals=(1,3))
    @test heisenberg"@(1,2)"==Heisenberg(centers=(1,2))
    @test heisenberg"ob(1,3)"==Heisenberg(orbitals=(1,3))
    @test heisenberg"sl(1,1)"==Heisenberg(atoms=(1,1))
    @test heisenberg""==Heisenberg()

    @test ising"x"==Ising('x') && ising"x sl(1,1)⊗ob(1,3)@(1,2)"==Ising('x',centers=(1,2),atoms=(1,1),orbitals=(1,3))
    @test ising"y"==Ising('y') && ising"y sl(1,1)⊗ob(1,3)@(1,2)"==Ising('y',centers=(1,2),atoms=(1,1),orbitals=(1,3))
    @test ising"z"==Ising('z') && ising"z sl(1,1)⊗ob(1,3)@(1,2)"==Ising('z',centers=(1,2),atoms=(1,1),orbitals=(1,3))

    @test gamma"x"==Gamma('x') && gamma"x sl(1,1)⊗ob(1,3)@(1,2)"==Gamma('x',centers=(1,2),atoms=(1,1),orbitals=(1,3))
    @test gamma"y"==Gamma('y') && gamma"y sl(1,1)⊗ob(1,3)@(1,2)"==Gamma('y',centers=(1,2),atoms=(1,1),orbitals=(1,3))
    @test gamma"z"==Gamma('z') && gamma"z sl(1,1)⊗ob(1,3)@(1,2)"==Gamma('z',centers=(1,2),atoms=(1,1),orbitals=(1,3))

    @test sˣ""==Sˣ() && sˣ"sl(1)⊗ob(2)"==Sˣ(atom=1,orbital=2)
    @test sʸ""==Sʸ() && sʸ"sl(1)⊗ob(2)"==Sʸ(atom=1,orbital=2)
    @test sᶻ""==Sᶻ() && sᶻ"sl(1)⊗ob(2)"==Sᶻ(atom=1,orbital=2)
end

@testset "couplings" begin
    @test @couplings(σ⁰"sp")==σ⁰("sp")
    @test @couplings(σ⁰"sp"+σᶻ"sp")==σ⁰("sp")+σᶻ("sp")
end
