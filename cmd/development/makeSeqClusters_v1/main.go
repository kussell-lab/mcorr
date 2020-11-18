package main

import (
	"fmt"
	"gopkg.in/alecthomas/kingpin.v2"
	"os"
	"os/exec"
	"path"
	"runtime"
	"time"
)

func main() {
	app := kingpin.New("makeSeqClusters_v1", "Takes the output of mcorr-pair, clusters sequences, and breaks clusters into core and flexible genomes")
	app.Version("v20201116")
	mcp := app.Arg("mcorr-pair_csv", "csv output from mcorr-pair as .csv file").Required().String()
	wrkdir := app.Arg("working_dir", "the working space and output directory").Required().String()
	msa := app.Arg("master_msa", "master msa file for all strain sequences").Required().String()
	cutoff := app.Flag("cutoff", "cutoff percentile of pairwise distances to use for making flat clusters (%)").Default("10").String()
	CFsplit := app.Flag("core-flex-split", "If you want to split genomes into core and flexible genes").Default("true").Bool()
	threshold := app.Flag("threshold", "threshold percentage above which you're considered a core gene (%)").Default("90").Int()
	ncpu := app.Flag("num-cpu", "Number of CPUs (default: using all available cores)").Default("0").Int()
	showProgress := app.Flag("show-progress", "Show progress").Default("false").Bool()
	kingpin.MustParse(app.Parse(os.Args[1:]))

	start := time.Now()

	//switch to the working directory
	os.Chdir(*wrkdir)
	//convert mcorr-pair output to distance matrix and list of strain names
	mcorrpairtodm(*mcp, *ncpu)

	//goExecPath, err := exec.LookPath("makeSeqClusters.py")
	//
	//if err != nil {
	//	fmt.Println("Error: ", err)
	//}else{
	//	fmt.Println("executable: ", goExecPath)
	//}
	////execute makeSeqClusters.py via commandline which will cluster the sequences
	////c := exec.Command("conda activate snakes")
	_, filename, _, ok := runtime.Caller(0)
	if !ok {
		panic("No caller information")
	}
	filepath := path.Join(path.Dir(filename), "makeSeqClusters.py")

	c := exec.Command(filepath, "distancematrix", "strains", *cutoff)
	if err := c.Run(); err != nil {
		fmt.Println("Error: ", err)
	}

	//cmd:= &exec.Cmd {
	//	Path: filepath,
	//	Args: []string{"python3", "makeSeqClusters.py", "distancematrix", "strains", *cutoff},
	//	Stdout: os.Stdout,
	//	Stderr: os.Stdout,
	//}
	//if err := cmd.Run(); err != nil {
	//	fmt.Println("Error:", err)
	//}

	//Step 3: make sequence cluster MSA files and core and flexible MSA files
	writeClusters(*msa, "cluster_list", *CFsplit, *threshold, *showProgress, *ncpu)

	duration := time.Since(start)
	fmt.Println(duration)
}
