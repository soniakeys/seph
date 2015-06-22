// Public domain.

package main

import (
	"bytes"
	"fmt"
	"io/ioutil"
	"log"
	"math"
	"strings"
	"time"

	"github.com/soniakeys/astro"
	"github.com/soniakeys/mpcformat"
	"github.com/soniakeys/observation"
	"github.com/soniakeys/sexagesimal"

	"github.com/soniakeys/meeus/base"
	"github.com/soniakeys/meeus/julian"
	"github.com/soniakeys/meeus/kepler"
	pp "github.com/soniakeys/meeus/planetposition"
	"github.com/soniakeys/meeus/solarxyz"

	"github.com/naoina/toml"
)

// part of all the data in the full orbit record.
type part struct {
	Desig, Epoch        string
	A, E                float64
	Inc, MA, Node, Peri float64
	G                   float64 `val:"defNaN"`
	H                   float64 `val:"defNaN"`
}

type job struct {
	Desig      string
	Start, End time.Time
	Site       string
}

func main() {
	// read job file
	td, err := ioutil.ReadFile("seph.job")
	if err != nil {
		log.Fatal(err)
	}
	var j job
	if err = toml.Unmarshal(td, &j); err != nil {
		log.Fatal(err)
	}
	// show stopper validations
	if j.Desig == "" {
		log.Fatal("no desig in seph.job")
	}
	if j.Start.IsZero() {
		log.Fatal("no start in seph.job")
	}
	var par *observation.ParallaxConst
	if strings.TrimSpace(j.Site) == "" {
		j.Site = "500"
	} else if j.Site != "500" {
		var p observation.ParallaxMap
		// bleh, parsing entire file wastes time here.
		p, err = mpcformat.ReadObscodeDatFile("obscode.dat")
		if err != nil {
			log.Fatal(err)
		}
		var ok bool
		if par, ok = p[j.Site]; !ok {
			log.Fatal("site ", j.Site, " unknown")
		}
	}
	_ = par
	// find desig in MPCORB.DAT
	bd, err := ioutil.ReadFile("MPCORB.DAT")
	if err != nil {
		log.Fatal(err)
	}
	lines := bytes.Split(bd, []byte{'\n'})
	var orbline []byte
	des := []byte(j.Desig)
	for _, line := range lines {
		if bytes.HasPrefix(line, des) {
			orbline = line
			break
		}
	}
	if len(orbline) == 0 {
		log.Fatal("desig ", j.Desig, " not found")
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
	// var inc time.Duration
	var el Elements
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
	sunPos := func(jde float64) (x, y, z, r float64) {
		// TODO replace with func that returns r as well)
		X, Y, Z := solarxyz.PositionJ2000(e, jde)
		return X, Y, Z, math.Sqrt(X*X + Y*Y + Z*Z)
	}

  fmt.Println()
	fmt.Println(j.Desig, "Epoch", y, m, d)
  fmt.Println("Site:", j.Site)
	fmt.Println("\n           Time      RA          Dec        V     Elongation")
	o := newOrbit(&el)
	α, δ, ψ, β, r, Δ := observation.AstrometricJ2000(julian.TimeToJD(j.Start),
		sunPos, o.Position)
	vs := ""
	if v := observation.Vmag(s.H, s.G, β, r, Δ); v >= 6 {
		vs = fmt.Sprintf("%4.1f", v)
	}
	printLine := func(t time.Time) {
		fmt.Printf("%s  %2v  %2v  %4s  %v\n", t.Format("2006-01-02 15:04:05"),
			sexa.NewFmtRA(α), sexa.NewFmtAngle(δ), vs, sexa.NewFmtAngle(ψ))
	}
	printLine(j.Start)
	if j.End.IsZero() {
		return
	}
	α, δ, ψ, β, r, Δ = observation.AstrometricJ2000(julian.TimeToJD(j.End),
		sunPos, o.Position)
	if vs > "" {
		vs = fmt.Sprintf("%4.1f", observation.Vmag(s.H, s.G, β, r, Δ))
	}
	printLine(j.End)
}

type Elements struct {
	Axis  float64 // Semimajor axis, a, in AU
	Ecc   float64 // Eccentricity, e
	Inc   float64 // Inclination, i, in radians
	ArgP  float64 // Argument of perihelion, ω, in radians
	Node  float64 // Longitude of ascending node, Ω, in radians
	TimeP float64 // Time of perihelion, T, as jde
}
type orbit struct {
	k          *Elements
	n          float64
	_A, _B, _C float64
	a, b, c    float64
}

func newOrbit(k *Elements) *orbit {
	o := &orbit{
		k: k,
		n: astro.K / k.Axis / math.Sqrt(k.Axis),
	}
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
	o._A = math.Atan2(F, P)
	o._B = math.Atan2(G, Q)
	o._C = math.Atan2(H, R)
	o.a = math.Hypot(F, P)
	o.b = math.Hypot(G, Q)
	o.c = math.Hypot(H, R)
	return o
}

func (o *orbit) Position(jde float64) (x, y, z, r float64) {
	M := o.n * (jde - o.k.TimeP)
	E, err := kepler.Kepler2b(o.k.Ecc, M, 15)
	if err != nil {
		E = kepler.Kepler3(o.k.Ecc, M)
	}
	r = kepler.Radius(E, o.k.Ecc, o.k.Axis)
	ν := kepler.True(E, o.k.Ecc)
	// (33.9) p. 229
	x = r * o.a * math.Sin(o._A+o.k.ArgP+ν)
	y = r * o.b * math.Sin(o._B+o.k.ArgP+ν)
	z = r * o.c * math.Sin(o._C+o.k.ArgP+ν)
	return
}
