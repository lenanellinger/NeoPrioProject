import os
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import gridspec
from helpers.get_data import get_feature_data, get_relevant_features
import sys
from optparse import OptionParser

output_dir = "/mnt/storage2/users/ahnelll1/master_thesis/analysis/images/cohorts"

def main(argv):
    usage = "usage: python cohorts_comparison_violin.py --majority-vote"
    desc = "Creates violin plot for pvacseq and neofox features."
    parser = OptionParser(usage=usage, description=desc)
    parser.add_option("--majority-vote", action="store_true", dest="majority_vote", default=False, help="If neoantigens should be filtered by majority vote of binding tools")
    (options, args) = parser.parse_args()

    majority_vote=options.majority_vote
    cohorts=['AxelMelanomaPhD', 'SomaticAndTreatment']
    
    features_data = get_feature_data(majority_vote, cohorts)

    relevant_features = get_relevant_features(True)

    fig = plt.figure(figsize=(30,8))
    gs = gridspec.GridSpec(1, len(relevant_features)) 

    for i, feature in enumerate(relevant_features):
        ax_violin = plt.subplot(gs[i])
        features_df = features_data.loc[:, [feature['name'], 'cohort']]
        features_df = features_df.rename(columns={feature['name']: 'value'})
        features_df['name'] = feature['name']
        sns.violinplot(ax=ax_violin, data=features_df, x="name", y='value', hue="cohort", split=True, inner="quart")
        if i != len(relevant_features) - 1:
            ax_violin.get_legend().set_visible(False)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'features_violin' + ('_majority_vote' if majority_vote else '') + '.png'))


if __name__ == "__main__":
    main(sys.argv)
