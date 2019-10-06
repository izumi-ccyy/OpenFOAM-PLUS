# createFields.H

define thermo variables from twoPhaseMixtureEThermo

# twoPhaseMixtureEThermo

## twoPhaseMixtureEThermo.H

```
Description
    Two phases thermo Internal energy mixture Defined as:
    e1 = Cv1(T - Tsat) + Hv1
    e2 = Cv2(T - Tsat) + Hv2
    e = (alpha1*rho1*e1 + alpha2*rho2*e2) / (alpha1*rho1 + alpha2*rho2)
```

$$
e_1 = C_{v1} ( T - T_{sat}) + h_{v1}
$$

$$
e_2 = C_{v2} ( T - T_{sat}) + h_{v2}
$$

$$
e = \frac{\alpha_1 \rho_1 e_1 + \alpha_2 \rho_2 e_2}{\alpha_1 \rho_1 + \alpha_2 \rho_2}
$$

## twoPhaseMixtureEThermo.C

```
void Foam::twoPhaseMixtureEThermo::correct() // correct T by alpha1 and alpha2
{
    incompressibleTwoPhaseMixture::correct();

    const volScalarField alpha1Rho1(alpha1()*rho1());
    const volScalarField alpha2Rho2(alpha2()*rho2());

    T_ =
        (
            (e_*(alpha1Rho1 + alpha2Rho2))
         -  (alpha1Rho1*Hf1() + alpha2Rho2*Hf2())
        )
       /(alpha1Rho1*Cv1() + alpha2Rho2*Cv2())
       + TSat_;

    T().correctBoundaryConditions();
}
```

$$
T = \frac{e (\alpha_1 \rho_1 + \alpha_2 \rho_2) - (\alpha_1 \rho_1 H_{f1} + \alpha_2 \rho_2 H_{f2})}{\alpha_1 \rho_1 C_{v1} + \alpha_2 \rho_2 C_{v2}} + T_{sat}
$$

where $H_f$ is the latent heat for phase

```
Foam::tmp<Foam::volScalarField> Foam::twoPhaseMixtureEThermo::he
(
    const volScalarField& p,
    const volScalarField& T
) const
{
    const volScalarField alpha1Rho1(alpha1()*rho1());
    const volScalarField alpha2Rho2(alpha2()*rho2());

    return
    (
        (T - TSat_)*(alpha1Rho1*Cv1() + alpha2Rho2*Cv2())
        + (alpha1Rho1*Hf1() + alpha2Rho2*Hf2())
    )
    / (alpha1Rho1 + alpha2Rho2);
}
```

$$
he = \frac{(T - T_{sat}) (\alpha_1 \rho_1 C_{v1} + \alpha_2 \rho_2 C_{v2}) }{\alpha_1 \rho_1 + \alpha_2 \rho_2}
$$

```
Foam::tmp<Foam::volScalarField> Foam::twoPhaseMixtureEThermo::hc() const // define hc
{
    const fvMesh& mesh = this->T_.mesh();

    return tmp<volScalarField>::New
    (
        IOobject
        (
            "hc",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("hc", Hf2() - Hf1())
    );
}
```

$$
h_c = H_{f2} - H_{f1}
$$

```
Foam::tmp<Foam::volScalarField> Foam::twoPhaseMixtureEThermo::Cp() const
{
    const volScalarField limitedAlpha1
    (
        min(max(alpha1_, scalar(0)), scalar(1))
    );

    return tmp<volScalarField>
    (
        new volScalarField
        (
            "cp",
            limitedAlpha1*Cp1() + (scalar(1) - limitedAlpha1)*Cp2()
        )
    );
}
```

$$
C_p = \alpha_1 C_{p1} + \alpha_2 C_{p2}
$$

where $\alpha_1$ and $\alpha_2$ are limited in $[0, 1]$

**This is the same for $\rho$, $C_v$ and $\kappa$** 

```
Foam::tmp<Foam::volScalarField> Foam::twoPhaseMixtureEThermo::gamma() const
{
    return tmp<volScalarField>
    (
        (alpha1_*Cp1() + alpha2_*Cp2())/(alpha1_*Cv1() + alpha2_*Cv2())
    );
}
```

$$
\gamma = \frac{\alpha_1 C_{p1} + \alpha_2 C_{p2}}{\alpha_1 C_{v1} + \alpha_2 C_{v2}}
$$
