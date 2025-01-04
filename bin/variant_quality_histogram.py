import pandas as pd
import matplotlib.pyplot as plt
import argparse


def parse_vcf(vcf_file):
    """Parse the VCF file and extract quality scores."""
    quality_scores = []
    with open(vcf_file, "r") as file:
        for line in file:
            if not line.startswith("#"):
                fields = line.strip().split("\t")
                quality_scores.append(float(fields[5]))
    return quality_scores


def plot_histogram(quality_scores, output_file):
    """Generate and save a histogram of variant quality scores."""
    plt.figure(figsize=(10, 6))
    plt.hist(quality_scores, bins=50, color="blue", edgecolor="black")
    plt.title("Histogram of Variant Quality Scores")
    plt.xlabel("Quality Score")
    plt.ylabel("Frequency")
    plt.grid(True)
    plt.savefig(output_file)
    plt.close()


def main():
    parser = argparse.ArgumentParser(
        description="Generate a histogram of variant quality scores from a VCF file."
    )
    parser.add_argument("vcf_file", help="Path to the input VCF file")
    parser.add_argument("output_file", help="Path to the output PNG file")
    args = parser.parse_args()

    quality_scores = parse_vcf(args.vcf_file)
    plot_histogram(quality_scores, args.output_file)


if __name__ == "__main__":
    main()
