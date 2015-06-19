// Public domain.

package main

import (
	"bytes"
	"fmt"
	"io/ioutil"
	"log"
	"math"
	"os"
	"time"

	"github.com/soniakeys/meeus/base"
	"github.com/soniakeys/meeus/elliptic"
	"github.com/soniakeys/meeus/julian"
	"github.com/soniakeys/meeus/kepler"
	pp "github.com/soniakeys/meeus/planetposition"
	"github.com/soniakeys/meeus/solarxyz"
	"github.com/soniakeys/mpcformat"
	"github.com/soniakeys/sexagesimal"
)

// part of all the data in the full orbit record.
type part struct {
	Desig, Epoch        string
	A, E                float64
	Inc, MA, Node, Peri float64
	G                   float64 `val:"defNaN"`
	H                   float64 `val:"defNaN"`
}


// astrometricJ2000 cut and paste from meeus/elliptic, modified to return
// additional values phase angle β and observer-object range Δ.
//
// AstrometricJ2000 is a utility function for computing astrometric coordinates.
//
// Argument f is a function that returns J2000 equatorial rectangular
// coodinates of a body.
//
// Results are J2000 right ascention, declination, and elongation.
func astrometricJ2000(f func(float64) (x, y, z, r float64), jde float64, e *pp.V87Planet) (α, δ, ψ, β, r, Δ float64) {
	X, Y, Z := solarxyz.PositionJ2000(e, jde)
	x, y, z, r := f(jde)
	// (33.10) p. 229
	ξ := X + x
	η := Y + y
	ζ := Z + z
	Δ = math.Sqrt(ξ*ξ + η*η + ζ*ζ)
	{
		τ := base.LightTime(Δ)
		x, y, z, r = f(jde - τ)
		ξ = X + x
		η = Y + y
		ζ = Z + z
		Δ = math.Sqrt(ξ*ξ + η*η + ζ*ζ)
	}
	α = math.Atan2(η, ξ)
	if α < 0 {
		α += 2 * math.Pi
	}
	δ = math.Asin(ζ / Δ)
	R := math.Sqrt(X*X + Y*Y + Z*Z)
	ψ = math.Acos((ξ*X + η*Y + ζ*Z) / R / Δ)
	β = math.Acos((ξ*x + η*y + ζ*z) / r / Δ)
	return
}

// position cut and paste from meeus/elliptic, modified to return additional
// values.
//
// Position returns observed equatorial coordinates of a body with Keplerian elements.
//
// Argument e must be a valid V87Planet object for Earth.
//
// Results are right ascension and declination α and δ, and elongation ψ,
// all in radians.  In addition, returns phase angle β, the sun-object distance
// r in AU and the earth-object distance Δ in AU.
func posFunc(k *elliptic.Elements, e *pp.V87Planet) func(jde float64) (α, δ, ψ, β, r, Δ float64) {
	// (33.6) p. 227
	n := base.K / k.Axis / math.Sqrt(k.Axis)
	const sε = base.SOblJ2000
	const cε = base.COblJ2000
	sΩ, cΩ := math.Sincos(k.Node)
	si, ci := math.Sincos(k.Inc)
	// (33.7) p. 228
	F := cΩ
	G := sΩ * cε
	H := sΩ * sε
	P := -sΩ * ci
	Q := cΩ*ci*cε - si*sε
	R := cΩ*ci*sε + si*cε
	// (33.8) p. 229
	A := math.Atan2(F, P)
	B := math.Atan2(G, Q)
	C := math.Atan2(H, R)
	a := math.Hypot(F, P)
	b := math.Hypot(G, Q)
	c := math.Hypot(H, R)

	f := func(jde float64) (x, y, z, r float64) {
		M := n * (jde - k.TimeP)
		E, err := kepler.Kepler2b(k.Ecc, M, 15)
		if err != nil {
			E = kepler.Kepler3(k.Ecc, M)
		}
		ν := kepler.True(E, k.Ecc)
		r = kepler.Radius(E, k.Ecc, k.Axis)
		// (33.9) p. 229
		x = r * a * math.Sin(A+k.ArgP+ν)
		y = r * b * math.Sin(B+k.ArgP+ν)
		z = r * c * math.Sin(C+k.ArgP+ν)
		return
	}
	return func(jde float64) (α, δ, ψ, β, r, Δ float64) {
		return astrometricJ2000(f, jde, e)
	}
}

func main() {
	des := []byte(os.Args[1])
	bd, err := ioutil.ReadFile("MPCORB.DAT")
	if err != nil {
		log.Fatal(err)
	}
	lines := bytes.Split(bd, []byte{'\n'})
	var orbline []byte
	for _, line := range lines {
		if bytes.HasPrefix(line, des) {
			orbline = line
			break
		}
	}
	if len(orbline) == 0 {
		log.Fatal(os.Args[1], " not found")
	}
	var s part
	su, err := mpcformat.NewExportUnmarshaler(&s)
	if err != nil {
		log.Fatal(err)
	}
	if err := su(orbline); err != nil {
		log.Fatal(err)
	}
	y, m, d, err := mpcformat.UnpackEpoch(s.Epoch)
	if err != nil {
		log.Fatal(err)
	}
  ep := julian.CalendarGregorianToJD(y, m, d)
	t1, err := time.Parse("2006-01-02T15:04:05Z", os.Args[2])
	if err != nil {
		log.Fatal(err)
	}
	t2 := t1
	// var inc time.Duration
	if len(os.Args) > 3 {
		t2, err = time.Parse("2006-01-02T15:04:05Z", os.Args[3])
		if err != nil {
			log.Fatal(err)
		}
	}
	var el elliptic.Elements
	el.Axis = s.A
	el.Ecc = s.E
	el.Inc = s.Inc * math.Pi / 180
	el.ArgP = s.Peri * math.Pi / 180
	el.Node = s.Node * math.Pi / 180
	el.TimeP = ep - (s.MA*math.Pi/180)*s.A*math.Sqrt(s.A)/base.K
	e, err := pp.LoadPlanet(pp.Earth)
	if err != nil {
		log.Fatal(err)
	}
	fmt.Println("  RA       Dec     V    Elongation")
	f := posFunc(&el, e)
	α, δ, ψ, β, r, Δ := f(julian.TimeToJD(t1))
	vs := ""
	if v := vmag(s.H, s.G, β, r, Δ); v >= 6 {
		vs = fmt.Sprintf("%4.1f", v)
	}
	printLine := func() {
		fmt.Printf("%v %v %4s, %v\n",
			sexa.NewFmtRA(α), sexa.NewFmtAngle(δ), vs, sexa.NewFmtAngle(ψ))
	}
	printLine()

	α, δ, ψ, β, r, Δ = f(julian.TimeToJD(t2))
	if vs > "" {
		vs = fmt.Sprintf("%4.1f", vmag(s.H, s.G, β, r, Δ))
	}
	printLine()
}

func vmag(H, G, β, r, Δ float64) float64 {
	if math.IsNaN(H) {
		return H
	}
	if math.IsNaN(G) {
		G = .15
	}
	tanβ2 := math.Tan(β / 2)
	Φ1 := math.Exp(-3.33 * math.Pow(tanβ2, .63))
	Φ2 := math.Exp(-1.87 * math.Pow(tanβ2, 1.22))
	return H + 5*math.Log10(r*Δ) - 2.5*math.Log10((1-G)*Φ1+G*Φ2)
}
