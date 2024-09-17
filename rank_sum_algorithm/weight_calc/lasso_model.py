import pandas as pd
import matplotlib.pyplot as plt

from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, roc_auc_score, roc_curve
from sklearn.inspection import permutation_importance
from sklearn import linear_model

from analysis.helpers.get_data import get_relevant_features_neofox
from data.get_data import get_train_data

class LassoModel:
    def __init__(self):
        self.train_features, self.train_labels, self.train_labels_int, self.test_features, self.test_labels, self.test_labels_int, self.test_patients = get_train_data()
        self.classifier = linear_model.Lasso(alpha=0.00001)
        self.classifier.fit(self.train_features, self.train_labels)

    def get_classifier(self):
        return self.classifier

    def save_feature_importances(self):
        feature_list = [name for name in get_relevant_features_neofox() if not name.startswith("Priority_score")]

        result = permutation_importance(
            self.classifier, self.test_features, self.test_labels_int, n_repeats=10, random_state=42, n_jobs=2
        )

        sorted_importances_idx = result.importances_mean.argsort()
        importances = pd.DataFrame(
            result.importances[sorted_importances_idx].T,
            columns=feature_list[sorted_importances_idx],
        )
        ax = importances.plot.box(vert=False, whis=10)
        ax.set_title("Permutation Importances (test set)")
        ax.axvline(x=0, color="k", linestyle="--")
        ax.set_xlabel("Decrease in accuracy score")
        ax.figure.tight_layout()

        plt.savefig(
            "/mnt/storage2/users/ahnelll1/master_thesis/NeoPrioProject/rank_sum_algorithm/weight_calc/images/lasso_feature_importances_test.png")


    def calc_statistics(self):
        pred_proba = self.classifier.predict(self.test_features)
        pred_labels = (pred_proba > 0.5).astype(int)

        precision = precision_score(self.test_labels_int, pred_labels)
        recall = recall_score(self.test_labels_int, pred_labels)
        f1 = f1_score(self.test_labels_int, pred_labels)
        auc = roc_auc_score(self.test_labels_int, pred_proba)

        print(f"LASSO train accuracy: {self.classifier.score(self.train_features, self.train_labels_int):.3f}")
        print(f"LASSO test accuracy: {self.classifier.score(self.test_features, self.test_labels_int):.3f}")
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
        plt.title('ROC Curve for LASSO')
        plt.legend()
        plt.savefig(
            "/mnt/storage2/users/ahnelll1/master_thesis/NeoPrioProject/rank_sum_algorithm/weight_calc/images/lasso_ROC.png")