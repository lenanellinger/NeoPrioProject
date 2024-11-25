import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc, rcParams
import seaborn as sns
import json
import sys

from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, roc_auc_score, roc_curve, precision_recall_curve, confusion_matrix, accuracy_score
from sklearn.inspection import permutation_importance
from sklearn.model_selection import cross_val_score, GridSearchCV, StratifiedKFold, train_test_split
from xgboost import XGBClassifier

from imblearn.under_sampling import RandomUnderSampler, TomekLinks
from imblearn.over_sampling import RandomOverSampler, SMOTE

sys.path.append(sys.path[0] + '/../../..')
from helpers.get_data import get_relevant_features_neofox
from data.get_data import get_train_data

rc('font', **{'family': 'serif', 'serif': ['cmr10'], 'size': 20})
rcParams['axes.unicode_minus'] = False


class XGBoostClassifierModel:
    def __init__(self):
        self.train_features_unchanged, self.train_labels_int_unchanged, self.train_features, self.train_labels_int, self.val_features, self.val_labels_int, self.test_features, self.test_labels_int, self.test_patients = get_train_data()
        
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
    
    def perform_cross_validation(self, sampler):
        skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
        
        fold_accuracies = []
        fold_recalls = []
        fold_precisions = []
        fold_f1s = []
        
        for train_index, val_index in skf.split(self.train_features_unchanged, self.train_labels_int_unchanged):
            X_train, X_eval = np.array(self.train_features_unchanged)[train_index], np.array(self.train_features_unchanged)[val_index]
            y_train, y_eval = np.array(self.train_labels_int_unchanged)[train_index], np.array(self.train_labels_int_unchanged)[val_index]
            
            X_eval_val, X_eval_test, y_eval_val, y_eval_test = train_test_split(X_eval, y_eval,
                                                                                test_size=0.5, random_state=42, stratify=y_eval)
            
            if sampler != None:
                X_train_res, y_train_res = sampler.fit_resample(X_train, y_train)
            else:
                X_train_res, y_train_res = X_train, y_train
            
            xgb = XGBClassifier(max_depth=3, scale_pos_weight=len([x for x in y_train_res if x == 0])/len([x for x in y_train_res if x == 1]), early_stopping_rounds=10, eval_metric="logloss", learning_rate=0.1, subsample=0.8, min_child_weight=20)

            xgb.fit(
                X_train_res, y_train_res,
                eval_set=[(X_eval_val, y_eval_val)],
                verbose=False
            )
            pred_eval = xgb.predict(X_eval_test)
            
            fold_accuracies.append(accuracy_score(pred_eval, y_eval_test))
            fold_recalls.append(recall_score(pred_eval, y_eval_test))
            fold_precisions.append(precision_score(pred_eval, y_eval_test))
            fold_f1s.append(f1_score(pred_eval, y_eval_test))

        print(f"Mean Accuracy: {np.mean(fold_accuracies):.4f} with std {np.std(fold_accuracies):.4f}")
        print(f"Mean Recall: {np.mean(fold_recalls):.4f} with std {np.std(fold_recalls):.4f}")
        print(f"Mean Precision: {np.mean(fold_precisions):.4f} with std {np.std(fold_precisions):.4f}")
        print(f"Mean F1 Score: {np.mean(fold_f1s):.4f} with std {np.std(fold_f1s):.4f}")
        
    def perform_over_under_sampling_model_selection(self):
        rus = RandomUnderSampler(sampling_strategy=0.5, random_state=42)
        tomek = TomekLinks()
        smote = SMOTE(sampling_strategy=0.5, random_state=42)
        ros = RandomOverSampler(sampling_strategy=0.5, random_state=42)
        
        for sampler, name in zip([None, rus, tomek, smote, ros], ['unchanged', 'RUS', 'Tomek', 'ROS', 'SMOTE']):
            print(name)
            self.perform_cross_validation(sampler)
    
    def get_classification_threshold(self):
        pred_proba = self.xgb_model.predict_proba(self.train_features)[:, 1]
        
        # ROC curve
        fpr, tpr, thresholds = roc_curve(self.train_labels_int, pred_proba)
        dif = tpr-fpr
        ind = dif.argmax()

        plt.plot(fpr, tpr, label='ROC curve (area = %0.2f)' % roc_auc_score(self.train_labels_int, pred_proba))
        plt.plot([0, 1], [0, 1], 'k--', label='No Skill')
        plt.plot(fpr[ind], tpr[ind], 'ro', label=f'best threshold {thresholds[ind]}')
        plt.xlim([0.0, 1.0])
        plt.ylim([0.0, 1.05])
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        plt.title('ROC Curve for XGBoost (train set)')
        plt.legend()
        plt.savefig("/mnt/storage2/users/ahnelll1/master_thesis/NeoPrioProject/prioritization/weight_calculation/images/xgb_ROC_train.png")
        
        # precision recall curve
        precision, recall, thresholds = precision_recall_curve(self.train_labels_int, pred_proba)

        # Plot Precision-Recall as a function of the threshold
        plt.plot(thresholds, precision[:-1], label="Precision")
        plt.plot(thresholds, recall[:-1], label="Recall")
        plt.xlabel("Threshold")
        plt.ylabel("Score")
        plt.title("Precision-Recall Trade-Off")
        plt.legend()
        plt.savefig("/mnt/storage2/users/ahnelll1/master_thesis/NeoPrioProject/prioritization/weight_calculation/images/xgb_precision_recall_train.png", dpi=300)
        
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
            self.xgb_model, self.test_features, self.test_labels_int, n_repeats=10, random_state=42, scoring='f1'
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
            if feature == 'Priority_score_imputed_fromRNA':
                weights[feature] = weights['Priority_score_imputed_fromDNA']
        total = 0
        for feature in get_relevant_features_neofox():
            if feature in weights:
                total += weights[feature]
            else:
                weights[feature] = np.mean(adjusted_importances)
                total += np.mean(adjusted_importances)

        for key in weights:
            weights[key] = weights[key] / total
        
        with open('/mnt/storage2/users/ahnelll1/master_thesis/NeoPrioProject/prioritization/data/weights.json', 'w') as f:
            json.dump(weights, f, ensure_ascii=False, indent=4)
            
        # weights horizontal bar chart
        plt.figure(figsize=(16, 1))
        sorted_weights = {k: v for k, v in sorted(weights.items(), key=lambda item: item[1])}
        weights_array = [weights[key] for key in sorted_weights]
        left = 0
        for label, w, c in zip(sorted_weights.keys(), weights_array, ["k", "#DD8047", "#A5AB81", "#D8B25C", "#854d2b", "#7BA79D", "#968C8C", "#94B6D2", "#b7bc9a", "#56756e", "#e4c98d"]):
            p = plt.barh("x", w, left=left, label=label, color=c)
            if w > 0.04:
                plt.bar_label(p, label_type='center', fmt='%.2f')
            elif w > 0.02:
                plt.bar_label(p, label_type='center', fmt='%.2f', rotation=90)                
            left += w
        plt.gca().axis('off')
        plt.savefig("/mnt/storage2/users/ahnelll1/master_thesis/NeoPrioProject/prioritization/weight_calculation/images/xgb_feature_importances_weights_prop_test.png", dpi=300)
        plt.close()
        
        sorted_importances_idx = result.importances_mean.argsort()
        importances = pd.DataFrame(
            result.importances[sorted_importances_idx].T,
            columns=[f.replace('_', ' ') for f in np.array(feature_list)[sorted_importances_idx]],
        )
        
        rc('font', **{'size': 20})
        ax = importances.plot.box(vert=False, whis=20, figsize=(12,6), colormap=sns.diverging_palette(166, 23, 45, 57, as_cmap=True))
        ax.axvline(x=0, color="k", linestyle="--")
        ax.set_xlabel("Decrease in F1 score")
        ax.figure.tight_layout()

        plt.savefig("/mnt/storage2/users/ahnelll1/master_thesis/NeoPrioProject/prioritization/weight_calculation/images/xgb_feature_importances_f1_scoring_test.png", dpi=300)
        plt.close()
        

    def calc_test_statistics(self):
        # Test Set
        pred_proba_test = self.xgb_model.predict_proba(self.test_features)[:, 1]
        pred_labels_test = (pred_proba_test > 0.6).astype(int)
        
        # Train Set
        pred_proba_train = self.xgb_model.predict_proba(self.train_features)[:, 1]
        pred_labels_train = (pred_proba_train > 0.6).astype(int)

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
        fpr_test, tpr_test, thresholds_test = roc_curve(self.test_labels_int, pred_proba_test)
        index_test = min(range(len(thresholds_test)), key=lambda i: abs(thresholds_test[i]- 0.6 ))
        
        fpr_train, tpr_train, thresholds_train = roc_curve(self.train_labels_int, pred_proba_train)
        index_train = min(range(len(thresholds_train)), key=lambda i: abs(thresholds_train[i]- 0.6 ))
        
        # Plot the ROC curve
        plt.figure(figsize=(8,6))
        plt.plot(fpr_test, tpr_test, label='Test set (area = %0.2f)' % roc_auc_score(self.test_labels_int, pred_proba_test), color="#94B6D2", linewidth=5)
        plt.plot(fpr_train, tpr_train, label='Train set (area = %0.2f)' % roc_auc_score(self.train_labels_int, pred_proba_train), color="#A5AB81", linewidth=5)
        plt.plot([0, 1], [0, 1], 'k--', label='No Skill')
        plt.plot(fpr_test[index_test], tpr_test[index_test], 'o', label=f'threshold 0.6', color="#DD8047", mew=5, ms=5)
        plt.plot(fpr_train[index_train], tpr_train[index_train], 'o', color="#DD8047", mew=5, ms=5)
        plt.xlim([0.0, 1.0])
        plt.ylim([0.0, 1.05])
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        plt.legend()
        plt.tight_layout()
        plt.savefig("/mnt/storage2/users/ahnelll1/master_thesis/NeoPrioProject/prioritization/weight_calculation/images/xgb_ROC.png", dpi=300)
        plt.close()
        
        cm = confusion_matrix(self.test_labels_int, pred_labels_test, labels=[0, 1])
        cm_df = pd.DataFrame(cm, columns=['N', 'P'], index=['N', 'P'])
        plt.figure(figsize=(6, 6))
        sns.heatmap(cm_df, annot=True, fmt='d', cmap=sns.light_palette("#7BA79D", as_cmap=True))
        plt.xlabel('Predicted Label')
        plt.ylabel('True Label')
        plt.savefig("/mnt/storage2/users/ahnelll1/master_thesis/NeoPrioProject/prioritization/weight_calculation/images/xgb_cm_test.png", dpi=300)
        plt.close()

        cm = confusion_matrix(self.train_labels_int, pred_labels_train, labels=[0, 1])
        plt.figure(figsize=(10, 7))
        sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', xticklabels=[0, 1], yticklabels=[0, 1])
        plt.xlabel('Predicted Label')
        plt.ylabel('True Label')
        plt.title('Confusion Matrix for XGBoost Classifier (train set)')
        plt.savefig("/mnt/storage2/users/ahnelll1/master_thesis/NeoPrioProject/prioritization/weight_calculation/images/xgb_cm_train.png")
        plt.close()

    def save_output(self):
        pred_proba = self.xgb_model.predict_proba(self.test_features)[:, 1]
        df = pd.DataFrame({'patientIdentifier': self.test_patients, 'pred_proba': pred_proba})
        df.to_csv('/mnt/storage2/users/ahnelll1/master_thesis/output_training_data/NEPdb_filtered_neofox_xgb_out.tsv', sep='\t', index=False, header=True)


if __name__ == '__main__':
    rf = XGBoostClassifierModel()
    rf.perform_parameter_grid_search()
    rf.perform_over_under_sampling_model_selection()
    rf.get_classification_threshold()
    rf.calc_test_statistics()
    rf.save_feature_importances()
    rf.save_output()