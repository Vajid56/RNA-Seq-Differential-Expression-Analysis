## ============================================================
##  Publication-standard Machine Learning for DEA data
##  Stroke vs Diabetes — Sayyad Vajid
##  Python equivalent of DEA_ML_Publication.R
## ============================================================

import pandas as pd
import numpy as np
import warnings
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
from scipy import stats
from scipy.stats import wilcoxon, kruskal, spearmanr, fisher_exact

from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.neighbors import KNeighborsClassifier
from sklearn.model_selection import (StratifiedKFold, cross_validate,
                                      train_test_split, GridSearchCV,
                                      learning_curve)
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.metrics import (confusion_matrix, ConfusionMatrixDisplay,
                              classification_report, roc_curve, auc,
                              precision_recall_curve, matthews_corrcoef,
                              f1_score, precision_score, recall_score,
                              accuracy_score)
from sklearn.pipeline import Pipeline

warnings.filterwarnings('ignore')
np.random.seed(42)

# ── Plot style ────────────────────────────────────────────────
plt.rcParams.update({
    'font.family'      : 'DejaVu Sans',
    'font.size'        : 11,
    'axes.spines.top'  : False,
    'axes.spines.right': False,
    'axes.labelsize'   : 12,
    'axes.titlesize'   : 13,
    'axes.titleweight' : 'bold',
    'figure.facecolor' : 'white',
    'axes.facecolor'   : '#FAFAFA',
    'grid.color'       : '#E0E0E0',
    'grid.linewidth'   : 0.7,
})

BLUE   = '#185FA5'
TEAL   = '#1D9E75'
CORAL  = '#D85A30'
AMBER  = '#BA7517'
PURPLE = '#534AB7'
LGRAY  = '#F5F5F3'

# ══════════════════════════════════════════════════════════════
# 1. LOAD DATA
# ══════════════════════════════════════════════════════════════
print("="*55)
print("  DEA Machine Learning Pipeline — Stroke vs Diabetes")
print("="*55)

df = pd.read_csv('DE_stroke_vs_diabetes.csv', index_col=0)
df.columns = ['baseMean', 'log2FoldChange', 'lfcSE',
              'stat', 'pvalue', 'padj']

print(f"\nLoaded: {df.shape[0]} genes, {df.shape[1]} columns")
print(df.head())

# ══════════════════════════════════════════════════════════════
# 2. CREATE CORRECT LABEL  (no data leakage)
# ══════════════════════════════════════════════════════════════
# Label = which direction the gene changes (biologically meaningful)
# NOT based on padj threshold (that causes 99.99% accuracy / leakage)

df['direction'] = np.where(df['log2FoldChange'] > 0,
                            'Stroke_higher', 'Diabetes_higher')

print("\nClass distribution:")
print(df['direction'].value_counts())
print(f"Class balance ratio: {df['direction'].value_counts().max() / df['direction'].value_counts().min():.2f}")

# ══════════════════════════════════════════════════════════════
# 3. FEATURE ENGINEERING  (safe — no leakage)
# ══════════════════════════════════════════════════════════════
# EXCLUDED: pvalue, padj (derived from stat — circular)
# INCLUDED: expression level, magnitude, uncertainty, signal/noise

df['log_baseMean']  = np.log1p(df['baseMean'])
df['abs_LFC']       = df['log2FoldChange'].abs()
df['lfc_se_ratio']  = df['abs_LFC'] / df['lfcSE']   # signal-to-noise
df['log_lfcSE']     = np.log1p(df['lfcSE'])

feature_cols = ['log_baseMean', 'abs_LFC', 'lfcSE',
                'log_lfcSE', 'stat', 'lfc_se_ratio']

X = df[feature_cols].values
y = df['direction'].values

print(f"\nFeatures used: {feature_cols}")
print(f"Missing values: {np.isnan(X).sum()}")

# ══════════════════════════════════════════════════════════════
# 4. TRAIN / TEST SPLIT  (stratified 80/20)
# ══════════════════════════════════════════════════════════════
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.2, random_state=42, stratify=y
)
print(f"\nTrain: {X_train.shape[0]} | Test: {X_test.shape[0]}")

# ══════════════════════════════════════════════════════════════
# 5. DEFINE MODELS WITH HYPERPARAMETER TUNING
# ══════════════════════════════════════════════════════════════
models = {
    'Random Forest': Pipeline([
        ('scaler', StandardScaler()),
        ('clf', GridSearchCV(
            RandomForestClassifier(random_state=42, n_jobs=-1),
            param_grid={'n_estimators': [100, 200],
                        'max_features': ['sqrt', 'log2'],
                        'max_depth'   : [None, 10, 20]},
            cv=3, scoring='f1', n_jobs=-1
        ))
    ]),
    'SVM (RBF)': Pipeline([
        ('scaler', StandardScaler()),
        ('clf', GridSearchCV(
            SVC(probability=True, random_state=42),
            param_grid={'C'    : [0.1, 1, 10],
                        'gamma': ['scale', 'auto', 0.1]},
            cv=3, scoring='f1', n_jobs=-1
        ))
    ]),
    'Logistic Regression': Pipeline([
        ('scaler', StandardScaler()),
        ('clf', GridSearchCV(
            LogisticRegression(max_iter=1000, random_state=42),
            param_grid={'C'      : [0.01, 0.1, 1, 10],
                        'penalty': ['l1', 'l2'],
                        'solver' : ['saga']},
            cv=3, scoring='f1', n_jobs=-1
        ))
    ]),
    'Gradient Boosting': Pipeline([
        ('scaler', StandardScaler()),
        ('clf', GridSearchCV(
            GradientBoostingClassifier(random_state=42),
            param_grid={'n_estimators'    : [100, 200],
                        'learning_rate'   : [0.05, 0.1],
                        'max_depth'       : [2, 3]},
            cv=3, scoring='f1', n_jobs=-1
        ))
    ]),
    'KNN': Pipeline([
        ('scaler', StandardScaler()),
        ('clf', GridSearchCV(
            KNeighborsClassifier(),
            param_grid={'n_neighbors': [3, 5, 7, 9, 11, 15]},
            cv=3, scoring='f1', n_jobs=-1
        ))
    ]),
}

# ══════════════════════════════════════════════════════════════
# 6. TRAIN ALL MODELS + 5-FOLD CROSS-VALIDATION
# ══════════════════════════════════════════════════════════════
cv5 = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

results      = {}
cv_results   = {}
trained_models = {}

print("\n" + "="*55)
print("  Training models (5-fold stratified CV)...")
print("="*55)

for name, pipeline in models.items():
    print(f"\n--- {name} ---")

    # Fit on full training set
    pipeline.fit(X_train, y_train)
    trained_models[name] = pipeline

    # Test set evaluation
    y_pred = pipeline.predict(X_test)
    y_prob = pipeline.predict_proba(X_test)[:, 1]

    # Encode for ROC
    le = LabelEncoder()
    y_test_enc = le.fit_transform(y_test)

    fpr, tpr, _ = roc_curve(y_test_enc, y_prob)
    roc_auc     = auc(fpr, tpr)
    mcc         = matthews_corrcoef(y_test, y_pred)

    results[name] = {
        'Accuracy'   : round(accuracy_score(y_test, y_pred) * 100, 2),
        'Precision'  : round(precision_score(y_test, y_pred,
                             pos_label='Stroke_higher') * 100, 2),
        'Recall'     : round(recall_score(y_test, y_pred,
                             pos_label='Stroke_higher') * 100, 2),
        'F1'         : round(f1_score(y_test, y_pred,
                             pos_label='Stroke_higher') * 100, 2),
        'AUC_ROC'    : round(roc_auc, 4),
        'MCC'        : round(mcc, 4),
        'fpr'        : fpr,
        'tpr'        : tpr,
    }

    # Cross-validation
    cv_sc = cross_validate(
        pipeline, X, y, cv=cv5,
        scoring=['accuracy', 'f1_weighted', 'roc_auc'],
        return_train_score=True, n_jobs=-1
    )
    cv_results[name] = cv_sc

    print(f"  Accuracy : {results[name]['Accuracy']}%")
    print(f"  F1       : {results[name]['F1']}%")
    print(f"  AUC-ROC  : {results[name]['AUC_ROC']}")
    print(f"  MCC      : {results[name]['MCC']}")

# ══════════════════════════════════════════════════════════════
# 7. RESULTS TABLE
# ══════════════════════════════════════════════════════════════
results_df = pd.DataFrame({
    k: {m: v for m, v in r.items()
        if m not in ('fpr', 'tpr')}
    for k, r in results.items()
}).T.reset_index()
results_df.columns = ['Model'] + list(results_df.columns[1:])

print("\n" + "="*55)
print("  FINAL RESULTS — HELD-OUT TEST SET")
print("="*55)
print(results_df.to_string(index=False))

results_df.to_csv('ML_results_publication.csv', index=False)
print("\nSaved: ML_results_publication.csv")

# ══════════════════════════════════════════════════════════════
# 8. PLOT 1 — ROC Curves
# ══════════════════════════════════════════════════════════════
colors_roc = [BLUE, TEAL, CORAL, AMBER, PURPLE]
fig, ax = plt.subplots(figsize=(8, 6))
ax.plot([0,1],[0,1],'k--', lw=1, alpha=0.4, label='Random baseline')

for (name, res), color in zip(results.items(), colors_roc):
    ax.plot(res['fpr'], res['tpr'], color=color, lw=2.2,
            label=f"{name}  (AUC = {res['AUC_ROC']:.4f})")
    ax.fill_between(res['fpr'], res['tpr'], alpha=0.04, color=color)

ax.set_xlabel('False Positive Rate')
ax.set_ylabel('True Positive Rate')
ax.set_title('ROC Curves — Stroke vs Diabetes Classification')
ax.legend(loc='lower right', frameon=False, fontsize=10)
ax.set_xlim([-0.02, 1.02]); ax.set_ylim([-0.02, 1.05])
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('roc_curves.png', dpi=300, bbox_inches='tight')
plt.close()
print("Saved: roc_curves.png")

# ══════════════════════════════════════════════════════════════
# 9. PLOT 2 — Confusion Matrices (all models)
# ══════════════════════════════════════════════════════════════
fig, axes = plt.subplots(2, 3, figsize=(15, 9))
fig.suptitle('Confusion Matrices — All Models (Held-out Test Set)',
             fontsize=14, fontweight='bold')
axes = axes.flatten()

for i, (name, pipeline) in enumerate(trained_models.items()):
    y_pred = pipeline.predict(X_test)
    cm     = confusion_matrix(y_test, y_pred,
                               labels=['Diabetes_higher', 'Stroke_higher'])
    disp   = ConfusionMatrixDisplay(
        cm, display_labels=['Diabetes\nhigher', 'Stroke\nhigher'])
    disp.plot(ax=axes[i], colorbar=False, cmap='Blues')
    axes[i].set_title(name, fontsize=11, fontweight='bold')
    for text in disp.text_.ravel():
        text.set_fontsize(14); text.set_fontweight('bold')

axes[-1].set_visible(False)
plt.tight_layout()
plt.savefig('confusion_matrices.png', dpi=300, bbox_inches='tight')
plt.close()
print("Saved: confusion_matrices.png")

# ══════════════════════════════════════════════════════════════
# 10. PLOT 3 — Feature Importance (Random Forest)
# ══════════════════════════════════════════════════════════════
rf_pipeline = trained_models['Random Forest']
rf_best     = rf_pipeline.named_steps['clf'].best_estimator_
importances = rf_best.feature_importances_

feat_df = pd.DataFrame({
    'Feature'   : feature_cols,
    'Importance': importances
}).sort_values('Importance', ascending=True)

feat_labels = {
    'log_baseMean' : 'log(baseMean)',
    'abs_LFC'      : '|log2FC|',
    'lfcSE'        : 'lfcSE',
    'log_lfcSE'    : 'log(lfcSE)',
    'stat'         : 'Wald statistic',
    'lfc_se_ratio' : '|log2FC| / lfcSE'
}
feat_df['Label'] = feat_df['Feature'].map(feat_labels)
feat_df['Color'] = feat_df['Importance'].apply(
    lambda v: BLUE if v > 0.2 else TEAL if v > 0.1 else AMBER)

fig, ax = plt.subplots(figsize=(9, 5))
bars = ax.barh(feat_df['Label'], feat_df['Importance'],
               color=feat_df['Color'], edgecolor='none', height=0.62)
for bar, val in zip(bars, feat_df['Importance']):
    ax.text(val + 0.003, bar.get_y() + bar.get_height()/2,
            f'{val*100:.1f}%', va='center', fontsize=10)
ax.set_xlabel('Feature Importance (Gini)')
ax.set_title('Random Forest Feature Importance\n'
             'Predicting expression direction: Stroke higher vs Diabetes higher')
ax.set_xlim(0, feat_df['Importance'].max() * 1.2)
ax.grid(axis='x', alpha=0.3)
plt.tight_layout()
plt.savefig('feature_importance.png', dpi=300, bbox_inches='tight')
plt.close()
print("Saved: feature_importance.png")

# ══════════════════════════════════════════════════════════════
# 11. PLOT 4 — CV scores comparison (boxplot per model)
# ══════════════════════════════════════════════════════════════
cv_auc_data = {name: cv_results[name]['test_roc_auc']
               for name in trained_models}
cv_df = pd.DataFrame(cv_auc_data)

fig, ax = plt.subplots(figsize=(9, 5))
bp = ax.boxplot(cv_df.values, patch_artist=True,
                medianprops=dict(color='white', linewidth=2))
for patch, color in zip(bp['boxes'], colors_roc):
    patch.set_facecolor(color); patch.set_alpha(0.8)
ax.set_xticklabels(cv_df.columns, rotation=15, ha='right')
ax.set_ylabel('AUC-ROC (5-fold CV)')
ax.set_title('Cross-Validation AUC-ROC — All Models')
ax.grid(axis='y', alpha=0.3)
plt.tight_layout()
plt.savefig('cv_comparison.png', dpi=300, bbox_inches='tight')
plt.close()
print("Saved: cv_comparison.png")

# ══════════════════════════════════════════════════════════════
# 12. PLOT 5 — Learning Curve (Random Forest)
# ══════════════════════════════════════════════════════════════
rf_simple = Pipeline([
    ('scaler', StandardScaler()),
    ('clf', RandomForestClassifier(n_estimators=100, random_state=42))
])
train_sizes, train_scores, test_scores = learning_curve(
    rf_simple, X, y, cv=5, scoring='f1_weighted',
    train_sizes=np.linspace(0.1, 1.0, 10), n_jobs=-1
)

fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(train_sizes, train_scores.mean(axis=1),
        color=BLUE, lw=2.2, marker='o', markersize=6, label='Train F1')
ax.plot(train_sizes, test_scores.mean(axis=1),
        color=CORAL, lw=2.2, marker='s', markersize=6, label='CV F1')
ax.fill_between(train_sizes,
                train_scores.mean(1) - train_scores.std(1),
                train_scores.mean(1) + train_scores.std(1),
                alpha=0.1, color=BLUE)
ax.fill_between(train_sizes,
                test_scores.mean(1) - test_scores.std(1),
                test_scores.mean(1) + test_scores.std(1),
                alpha=0.1, color=CORAL)
ax.set_xlabel('Training set size')
ax.set_ylabel('F1 Score (weighted)')
ax.set_title('Learning Curve — Random Forest')
ax.legend(frameon=False); ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('learning_curve.png', dpi=300, bbox_inches='tight')
plt.close()
print("Saved: learning_curve.png")

# ══════════════════════════════════════════════════════════════
# 13. STATISTICAL TESTS FOR PUBLICATION
# ══════════════════════════════════════════════════════════════
print("\n" + "="*55)
print("  STATISTICAL TESTS")
print("="*55)

sig   = df[(df['padj'] < 0.05) & (df['log2FoldChange'].abs() > 1)]
up    = sig[sig['log2FoldChange'] >  1]
down  = sig[sig['log2FoldChange'] < -1]

# Test 1: Wilcoxon rank-sum (Mann-Whitney U)
print("\n--- Wilcoxon rank-sum: baseMean (Up vs Down genes) ---")
stat_w, p_w = stats.mannwhitneyu(
    up['baseMean'], down['baseMean'], alternative='two-sided')
r_effect = stat_w / (len(up) * len(down))
print(f"  Statistic : {stat_w:.2f}")
print(f"  p-value   : {p_w:.4e}")
print(f"  Effect (r): {r_effect:.4f}")

# Test 2: Kruskal-Wallis across significance categories
print("\n--- Kruskal-Wallis: |LFC| across significance categories ---")
df['category'] = 'Not_sig'
df.loc[(df['padj'] < 0.001) & (df['log2FoldChange'].abs() > 1),
       'category'] = 'Highly_sig'
df.loc[(df['padj'] < 0.05)  & (df['log2FoldChange'].abs() > 1) &
       (df['category'] == 'Not_sig'), 'category'] = 'Significant'

groups = [g['log2FoldChange'].abs().values
          for _, g in df.groupby('category')]
stat_k, p_k = kruskal(*groups)
print(f"  H-statistic: {stat_k:.4f}")
print(f"  p-value    : {p_k:.4e}")

# Post-hoc: Dunn test (manual pairwise Mann-Whitney + Bonferroni)
print("\n--- Pairwise Mann-Whitney (Bonferroni corrected) ---")
cats = df['category'].unique()
n_comparisons = len(cats) * (len(cats)-1) // 2
for i, c1 in enumerate(cats):
    for c2 in cats[i+1:]:
        g1 = df[df['category']==c1]['log2FoldChange'].abs()
        g2 = df[df['category']==c2]['log2FoldChange'].abs()
        _, p = stats.mannwhitneyu(g1, g2, alternative='two-sided')
        p_adj = min(p * n_comparisons, 1.0)
        sig_star = '***' if p_adj<0.001 else '**' if p_adj<0.01 else '*' if p_adj<0.05 else 'ns'
        print(f"  {c1} vs {c2}: p_adj={p_adj:.4e} {sig_star}")

# Test 3: Spearman correlation
print("\n--- Spearman correlation: log(baseMean) vs |LFC| ---")
rho, p_sp = spearmanr(np.log1p(df['baseMean']),
                       df['log2FoldChange'].abs())
print(f"  rho     : {rho:.4f}")
print(f"  p-value : {p_sp:.4e}")

# Test 4: Fisher's exact test
print("\n--- Fisher's exact test: high expression vs upregulation ---")
high_expr   = df['baseMean'] > df['baseMean'].median()
upregulated = (df['log2FoldChange'] > 1) & (df['padj'] < 0.05)
ct = pd.crosstab(high_expr, upregulated)
odds, p_fi = fisher_exact(ct)
print(f"  Odds ratio: {odds:.4f}")
print(f"  p-value   : {p_fi:.4e}")

# ══════════════════════════════════════════════════════════════
# 14. PLOT 6 — Full summary dashboard
# ══════════════════════════════════════════════════════════════
fig = plt.figure(figsize=(16, 12))
fig.suptitle('ML Analysis Dashboard — Stroke vs Diabetes DEA\n'
             'Sayyad Vajid | University of Hyderabad',
             fontsize=15, fontweight='bold', y=0.98)
gs = gridspec.GridSpec(2, 3, figure=fig, hspace=0.42, wspace=0.35)

# ROC
ax1 = fig.add_subplot(gs[0, 0])
ax1.plot([0,1],[0,1],'k--', lw=0.8, alpha=0.4)
for (name, res), color in zip(results.items(), colors_roc):
    ax1.plot(res['fpr'], res['tpr'], color=color, lw=1.8,
             label=f"{name.split()[0]} ({res['AUC_ROC']:.3f})")
ax1.set_title('ROC Curves'); ax1.set_xlabel('FPR'); ax1.set_ylabel('TPR')
ax1.legend(frameon=False, fontsize=8)
ax1.set_xlim([-0.02,1.02]); ax1.set_ylim([-0.02,1.05])
ax1.grid(True, alpha=0.25)

# Best model confusion matrix (RF)
ax2 = fig.add_subplot(gs[0, 1])
y_pred_rf = trained_models['Random Forest'].predict(X_test)
cm_rf = confusion_matrix(y_test, y_pred_rf,
                          labels=['Diabetes_higher','Stroke_higher'])
sns.heatmap(cm_rf, annot=True, fmt='d', cmap='Blues', ax=ax2,
            xticklabels=['Diabetes\nhigher','Stroke\nhigher'],
            yticklabels=['Diabetes\nhigher','Stroke\nhigher'],
            annot_kws={'size':13,'weight':'bold'}, cbar=False, linewidths=1)
ax2.set_title('Confusion Matrix\n(Random Forest)')
ax2.set_xlabel('Predicted'); ax2.set_ylabel('True')

# Feature importance
ax3 = fig.add_subplot(gs[0, 2])
feat_sorted = feat_df.sort_values('Importance', ascending=True)
ax3.barh(feat_sorted['Label'], feat_sorted['Importance'],
         color=feat_sorted['Color'], edgecolor='none', height=0.6)
ax3.set_title('Feature Importance\n(Random Forest)')
ax3.set_xlabel('Importance (Gini)'); ax3.grid(axis='x', alpha=0.3)

# CV boxplot
ax4 = fig.add_subplot(gs[1, 0])
bp2 = ax4.boxplot(cv_df.values, patch_artist=True,
                   medianprops=dict(color='white', lw=2))
for patch, color in zip(bp2['boxes'], colors_roc):
    patch.set_facecolor(color); patch.set_alpha(0.8)
ax4.set_xticklabels([n.split()[0] for n in cv_df.columns],
                     rotation=20, ha='right', fontsize=9)
ax4.set_ylabel('AUC-ROC'); ax4.set_title('5-Fold CV AUC-ROC')
ax4.grid(axis='y', alpha=0.3)

# Learning curve
ax5 = fig.add_subplot(gs[1, 1])
ax5.plot(train_sizes, train_scores.mean(1), color=BLUE, lw=2,
         marker='o', markersize=5, label='Train')
ax5.plot(train_sizes, test_scores.mean(1), color=CORAL, lw=2,
         marker='s', markersize=5, label='CV')
ax5.fill_between(train_sizes, train_scores.mean(1)-train_scores.std(1),
                 train_scores.mean(1)+train_scores.std(1), alpha=0.1, color=BLUE)
ax5.fill_between(train_sizes, test_scores.mean(1)-test_scores.std(1),
                 test_scores.mean(1)+test_scores.std(1), alpha=0.1, color=CORAL)
ax5.set_title('Learning Curve (RF)'); ax5.set_xlabel('Train size')
ax5.set_ylabel('F1 Score'); ax5.legend(frameon=False, fontsize=9)
ax5.grid(True, alpha=0.25)

# Metrics bar chart
ax6 = fig.add_subplot(gs[1, 2])
metrics = ['Accuracy', 'F1', 'AUC_ROC']
x = np.arange(len(metrics))
w = 0.15
for i, (name, color) in enumerate(zip(results.keys(), colors_roc)):
    vals = [results[name]['Accuracy']/100,
            results[name]['F1']/100,
            results[name]['AUC_ROC']]
    ax6.bar(x + i*w, vals, w, label=name.split()[0],
            color=color, alpha=0.85, edgecolor='none')
ax6.set_xticks(x + w*2)
ax6.set_xticklabels(['Accuracy', 'F1', 'AUC-ROC'])
ax6.set_ylim([0.5, 1.05]); ax6.set_title('Metrics Comparison')
ax6.legend(frameon=False, fontsize=8, ncol=2)
ax6.grid(axis='y', alpha=0.3)
ax6.yaxis.set_major_formatter(plt.FuncFormatter(lambda v,_: f'{v:.2f}'))

plt.savefig('ML_dashboard.png', dpi=300, bbox_inches='tight')
plt.close()
print("\nSaved: ML_dashboard.png")

# ══════════════════════════════════════════════════════════════
# 15. SAVE BEST MODEL
# ══════════════════════════════════════════════════════════════
import joblib
best_name = results_df.loc[results_df['AUC_ROC'].astype(float).idxmax(), 'Model']
joblib.dump(trained_models[best_name], 'best_model.pkl')
print(f"\nBest model: {best_name}")
print("Saved: best_model.pkl")

print("\n" + "="*55)
print("  ALL OUTPUTS SAVED")
print("="*55)
print("  ML_results_publication.csv")
print("  roc_curves.png")
print("  confusion_matrices.png")
print("  feature_importance.png")
print("  cv_comparison.png")
print("  learning_curve.png")
print("  ML_dashboard.png")
print("  best_model.pkl")
print("\n  Done!")
