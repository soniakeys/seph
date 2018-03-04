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
	"github.com/soniakeys/unit"

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
	ep := astro.FFCalendarGregorianToJD(y, m, d)
	// var inc time.Duration
	var el astro.Elements
	el.Axis = s.A
	el.Ecc = s.E
	el.Inc = unit.AngleFromDeg(s.Inc)
	el.ArgP = unit.AngleFromDeg(s.Peri)
	el.Node = unit.AngleFromDeg(s.Node)
	el.TimeP = ep - (s.MA*math.Pi/180)*s.A*math.Sqrt(s.A)/astro.K

	e, err := astro.LoadPlanet(astro.Earth)
	if err != nil {
		log.Fatal(err)
	}
	sunPos := func(jde float64) (x, y, z, r float64) {
		return astro.SolarPositionJ2000(e, jde)
	}

	fmt.Println()
	fmt.Println(s.Desig, "Epoch", y, m, d)
	fmt.Println("Site:", j.Site)
	fmt.Println("\n           Time      RA          Dec        V     Elongation")
	o := astro.NewOrbit(&el)
	α, δ, ψ, β, r, Δ := observation.AstrometricJ2000(astro.TimeToJD(j.Start),
		sunPos, o.Position)
	vs := ""
	if v := observation.Vmag(s.H, s.G, β, r, Δ); v >= 6 {
		vs = fmt.Sprintf("%4.1f", v)
	}
	printLine := func(t time.Time) {
		fmt.Printf("%s  %2v  %2v  %4s  %v\n", t.Format("2006-01-02 15:04:05"),
			sexa.FmtRA(unit.RA(α)),
			sexa.FmtAngle(unit.Angle(δ)),
			vs,
			sexa.FmtAngle(unit.Angle(ψ)))
	}
	printLine(j.Start)
	if j.End.IsZero() {
		return
	}
	α, δ, ψ, β, r, Δ = observation.AstrometricJ2000(astro.TimeToJD(j.End),
		sunPos, o.Position)
	if vs > "" {
		vs = fmt.Sprintf("%4.1f", observation.Vmag(s.H, s.G, β, r, Δ))
	}
	printLine(j.End)
}
