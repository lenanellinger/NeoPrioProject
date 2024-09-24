import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import json
import sys

from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, roc_auc_score, roc_curve
from sklearn.inspection import permutation_importance

sys.path.append(sys.path[0] + '/../../analysis')
from helpers.get_data import get_relevant_features_neofox
sys.path.append(sys.path[0] + '/../../rank_sum_algorithm')
from weight_calc.data.get_data import get_train_data


class RandomForestModel:
    def __init__(self):
        self.train_features, self.train_labels, self.train_labels_int, self.test_features, self.test_labels, self.test_labels_int, self.test_patients = get_train_data()
        self.classifier = RandomForestClassifier(n_estimators=100, random_state=42, min_samples_leaf=20)
        self.classifier.fit(self.train_features, self.train_labels)

    def get_classifier(self):
        return self.classifier

    def save_feature_importances(self):
        feature_list = [name for name in get_relevant_features_neofox() if not name.startswith("Priority_score")]

        result = permutation_importance(
            self.classifier, self.test_features, self.test_labels, n_repeats=10, random_state=42, n_jobs=2
        )

        print("Most Important Features:")
        for i in result.importances_mean.argsort()[::-1]:
            if result.importances_mean[i] - 2 * result.importances_std[i] > 0:
                print(f"\t{feature_list[i]:<8}"
                      f"{result.importances_mean[i]:.3f}"
                      f" +/- {result.importances_std[i]:.3f}")

        importances_mean = result.importances_mean
        importances_std = result.importances_std

        # Adjust the importance scores by dividing mean by std (to penalize high variability)
        adjusted_importances = importances_mean / (importances_std + 1e-10)  # Adding a small value to avoid division by zero

        # Normalize the adjusted importances to get final weights
        feature_weights_adjusted = adjusted_importances / adjusted_importances.sum()

        weights = {}
        for feature, weight in zip(feature_list, feature_weights_adjusted):
            weights[feature] = weight if weight > 0 else 1e-10
        total = 1
        for feature in get_relevant_features_neofox():
            if feature not in weights:
                weights[feature] = np.mean(feature_weights_adjusted)
                total += np.mean(feature_weights_adjusted)

        for key in weights:
            weights[key] = weights[key] / total

        with open('/mnt/storage2/users/ahnelll1/master_thesis/NeoPrioProject/rank_sum_algorithm/data/weights.json', 'w') as f:
            json.dump(weights, f, ensure_ascii=False, indent=4)

        sorted_importances_idx = result.importances_mean.argsort()
        importances = pd.DataFrame(
            result.importances[sorted_importances_idx].T,
            columns=np.array(feature_list)[sorted_importances_idx],
        )
        ax = importances.plot.box(vert=False, whis=10)
        ax.set_title("Permutation Importances (test set)")
        ax.axvline(x=0, color="k", linestyle="--")
        ax.set_xlabel("Decrease in accuracy score")
        ax.figure.tight_layout()

        plt.savefig("/mnt/storage2/users/ahnelll1/master_thesis/NeoPrioProject/rank_sum_algorithm/weight_calc/images/rf_feature_importances_test.png")

    def calc_test_statistics(self):
        # Test Set
        pred_proba = self.classifier.predict_proba(self.test_features)[:, 1]
        pred_labels = (pred_proba > 0.5).astype(int)

        precision = precision_score(self.test_labels_int, pred_labels)
        recall = recall_score(self.test_labels_int, pred_labels)
        f1 = f1_score(self.test_labels_int, pred_labels)
        auc = roc_auc_score(self.test_labels_int, pred_proba)

        print(f"RF train accuracy: {self.classifier.score(self.train_features, self.train_labels):.3f}")
        print(f"RF test accuracy: {self.classifier.score(self.test_features, self.test_labels):.3f}")
        print("Precision:", precision)
        print("Recall:", recall)
        print("F1-score:", f1)
        print("AUC:", auc)

        # Calculate ROC curve
        fpr, tpr, thresholds = roc_curve(self.test_labels_int, pred_proba)
        # Plot the ROC curve
        plt.figure()
        plt.plot(fpr, tpr, label='ROC curve (area = %0.2f)' % auc)
        plt.plot([0, 1], [0, 1], 'k--', label='No Skill')
        plt.xlim([0.0, 1.0])
        plt.ylim([0.0, 1.05])
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        plt.title('ROC Curve for Random Forest')
        plt.legend()
        plt.savefig("/mnt/storage2/users/ahnelll1/master_thesis/NeoPrioProject/rank_sum_algorithm/weight_calc/images/rf_ROC.png")


    def save_output(self):
        pred_proba = self.classifier.predict_proba(self.test_features)[:, 1]
        df = pd.DataFrame({'patientIdentifier': self.test_patients, 'pred_proba': pred_proba})
        df.to_csv('/mnt/storage2/users/ahnelll1/master_thesis/output_training_data/NEPdb_neofox_annotations_random_forest_out.tsv', sep='\t', index=False, header=True)


if __name__ == '__main__':
    rf = RandomForestModel()
    rf.save_feature_importances()
    rf.calc_test_statistics()
    rf.save_output()