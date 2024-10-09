<?php
ini_set("display_errors",1);
error_reporting(E_ALL);
$netid=$_POST["netid"];
$dataset=$_POST["name"];
$commands=$_POST["commands"];
$fname="Madgraph-Pythia-Delphes-Data/$netid-$dataset.inp";
$file = fopen($fname, "w") or die("Unable to open file $fname!");
fwrite($file, "$commands");
echo "Thank you";
?>
