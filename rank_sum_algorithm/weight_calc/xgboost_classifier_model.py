import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import json
import sys

from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, roc_auc_score, roc_curve, precision_recall_curve, confusion_matrix, accuracy_score
from sklearn.inspection import permutation_importance
from sklearn.model_selection import cross_val_score, GridSearchCV, StratifiedKFold
from xgboost import XGBClassifier

sys.path.append(sys.path[0] + '/../../analysis')
from helpers.get_data import get_relevant_features_neofox
sys.path.append(sys.path[0] + '/../../rank_sum_algorithm')
from weight_calc.data.get_data import get_train_data


class XGBoostClassifierModel:
    def __init__(self):
        self.train_features, self.train_labels_int, self.val_features, self.val_labels_int, self.test_features, self.test_labels_int, self.test_patients = get_train_data()
        
        eval_set = [(self.val_features, self.val_labels_int)]
        self.xgb_model = XGBClassifier(max_depth=3, scale_pos_weight=len([x for x in self.train_labels_int if x == 0])/len([x for x in self.train_labels_int if x == 1]), early_stopping_rounds=10, eval_metric="logloss", learning_rate=0.1, subsample=0.8, min_child_weight=20)
        self.xgb_model.fit(self.train_features, self.train_labels_int, eval_set=eval_set, verbose=True)
        
    def perform_parameter_grid_search(self):
        param_grid = {
            'min_child_weight': [1, 10, 20],
            'subsample': [0.3, 0.5, 0.8],
            'max_depth': [3, 4, 5],
        }
        
        xgb = XGBClassifier(early_stopping_rounds=10, eval_metric="logloss", learning_rate=0.1, scale_pos_weight=len([x for x in self.train_labels_int if x == 0])/len([x for x in self.train_labels_int if x == 1]))
                    
        skf = StratifiedKFold(n_splits=3, shuffle=True, random_state=42)
        grid_search = GridSearchCV(xgb, param_grid=param_grid, scoring='f1', cv=skf.split(self.train_features, self.train_labels_int))

        grid_search.fit(self.train_features, self.train_labels_int, eval_set=[(self.val_features, self.val_labels_int)])

        print("Best Parameters:", grid_search.best_params_)        
    
    def perform_cross_validation(self):
        skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
        
        fold_accuracies = []
        
        for train_index, val_index in skf.split(self.train_features, self.train_labels_int):
            X_train, X_eval = np.array(self.train_features)[train_index], np.array(self.train_features)[val_index]
            y_train, y_eval = np.array(self.train_labels_int)[train_index], np.array(self.train_labels_int)[val_index]
            
            xgb = XGBClassifier(max_depth=3, scale_pos_weight=len([x for x in self.train_labels_int if x == 0])/len([x for x in self.train_labels_int if x == 1]), early_stopping_rounds=10, eval_metric="logloss", learning_rate=0.1, subsample=0.8, min_child_weight=20)

            xgb.fit(
                X_train, y_train,
                eval_set=[(self.val_features, self.val_labels_int)]
            )
            
            accuracy = xgb.score(X_eval, y_eval)
            fold_accuracies.append(accuracy)

        print(f"Mean Accuracy: {np.mean(fold_accuracies):.4f}")
        print(f"Standard Deviation: {np.std(fold_accuracies):.4f}")
    
    def get_classification_threshold(self):
        pred_proba = self.xgb_model.predict_proba(self.train_features)[:, 1]
        
        # precision recall curve
        precision, recall, thresholds = precision_recall_curve(self.train_labels_int, pred_proba)

        # Plot Precision-Recall as a function of the threshold
        plt.plot(thresholds, precision[:-1], label="Precision")
        plt.plot(thresholds, recall[:-1], label="Recall")
        plt.xlabel("Threshold")
        plt.ylabel("Score")
        plt.title("Precision-Recall Trade-Off")
        plt.legend()
        plt.savefig("/mnt/storage2/users/ahnelll1/master_thesis/NeoPrioProject/rank_sum_algorithm/weight_calc/images/xgb_precision_recall_train.png", dpi=300)
        
        # best F1 score
        best_f1 = 0
        best_threshold = 0

        for threshold in np.arange(0, 1, 0.1).tolist():
            y_pred = (pred_proba >= threshold).astype(int)
            f1 = f1_score(self.train_labels_int, y_pred)
            if f1 > best_f1:
                best_f1 = f1
                best_threshold = threshold

        print("Best Threshold:", best_threshold)
        print("Best F1 Score:", best_f1)
        
        # Youden’s J Statistic
        fpr, tpr, thresholds = roc_curve(self.train_labels_int, pred_proba)

        youdens_j = tpr - fpr
        best_threshold = thresholds[youdens_j.argmax()]

        print("Best Threshold based on Youden’s J:", best_threshold)
        

    def save_feature_importances(self):
        feature_list = get_relevant_features_neofox()

        result = permutation_importance(
            self.xgb_model, self.test_features, self.test_labels_int, n_repeats=10, random_state=42
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

        weights = {}
        for feature, weight in zip(feature_list, adjusted_importances):
            weights[feature] = weight if weight > 0 else 1e-10
        total = 0
        for feature in get_relevant_features_neofox():
            if feature in weights:
                total += weights[feature]
            else:
                weights[feature] = np.mean(adjusted_importances)
                total += np.mean(adjusted_importances)

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

        plt.savefig("/mnt/storage2/users/ahnelll1/master_thesis/NeoPrioProject/rank_sum_algorithm/weight_calc/images/xgb_feature_importances_test.png")
        

    def calc_test_statistics(self):
        # Test Set
        pred_proba_test = self.xgb_model.predict_proba(self.test_features)[:, 1]
        pred_labels_test = (pred_proba_test > 0.7).astype(int)
        
        # Train Set
        pred_proba_train = self.xgb_model.predict_proba(self.train_features)[:, 1]
        pred_labels_train = (pred_proba_train > 0.7).astype(int)

        print("XGB test accuracy:", accuracy_score(self.test_labels_int, pred_labels_test))
        print("XGB test Precision:", precision_score(self.test_labels_int, pred_labels_test))
        print("XGB test Recall:", recall_score(self.test_labels_int, pred_labels_test))
        print("XGB test F1-score:", f1_score(self.test_labels_int, pred_labels_test))
        print("XGB test AUC:", roc_auc_score(self.test_labels_int, pred_proba_test))
        print("\nXGB train accuracy:", accuracy_score(self.train_labels_int, pred_labels_train))
        print("XGB train Precision:", precision_score(self.train_labels_int, pred_labels_train))
        print("XGB train Recall:", recall_score(self.train_labels_int, pred_labels_train))
        print("XGB train F1-score:", f1_score(self.train_labels_int, pred_labels_train))
        print("XGB train AUC:", roc_auc_score(self.train_labels_int, pred_proba_train))

        # Calculate ROC curve
        fpr, tpr, thresholds = roc_curve(self.test_labels_int, pred_proba_test)
        # Plot the ROC curve
        plt.figure()
        plt.plot(fpr, tpr, label='ROC curve (area = %0.2f)' % roc_auc_score(self.test_labels_int, pred_proba_test))
        plt.plot([0, 1], [0, 1], 'k--', label='No Skill')
        plt.xlim([0.0, 1.0])
        plt.ylim([0.0, 1.05])
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        plt.title('ROC Curve for XGBoost (test set)')
        plt.legend()
        plt.savefig("/mnt/storage2/users/ahnelll1/master_thesis/NeoPrioProject/rank_sum_algorithm/weight_calc/images/xgb_ROC.png")
        
        cm = confusion_matrix(self.test_labels_int, pred_labels_test, labels=[0, 1])
        plt.figure(figsize=(10, 7))
        sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', xticklabels=[0, 1], yticklabels=[0, 1])
        plt.xlabel('Predicted Label')
        plt.ylabel('True Label')
        plt.title('Confusion Matrix for XGBoost Classifier (test set)')
        plt.savefig("/mnt/storage2/users/ahnelll1/master_thesis/NeoPrioProject/rank_sum_algorithm/weight_calc/images/xgb_cm_test.png")

        cm = confusion_matrix(self.train_labels_int, pred_labels_train, labels=[0, 1])
        plt.figure(figsize=(10, 7))
        sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', xticklabels=[0, 1], yticklabels=[0, 1])
        plt.xlabel('Predicted Label')
        plt.ylabel('True Label')
        plt.title('Confusion Matrix for XGBoost Classifier (train set)')
        plt.savefig("/mnt/storage2/users/ahnelll1/master_thesis/NeoPrioProject/rank_sum_algorithm/weight_calc/images/xgb_cm_train.png")

    def save_output(self):
        pred_proba = self.xgb_model.predict_proba(self.test_features)[:, 1]
        df = pd.DataFrame({'patientIdentifier': self.test_patients, 'pred_proba': pred_proba})
        df.to_csv('/mnt/storage2/users/ahnelll1/master_thesis/output_training_data/NEPdb_neofox_annotations_xgb_out.tsv', sep='\t', index=False, header=True)


if __name__ == '__main__':
    rf = XGBoostClassifierModel()
    #rf.perform_parameter_grid_search()
    #rf.perform_cross_validation()
    #rf.get_classification_threshold()
    #rf.save_feature_importances()
    rf.calc_test_statistics()
    #rf.save_output()