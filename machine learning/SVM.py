import pandas as pd
import numpy as np
from sklearn.svm import SVC
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score, precision_recall_curve, auc, RocCurveDisplay, PrecisionRecallDisplay
import matplotlib.pyplot as plt

# 读取数据
file_path = "C:/Users/31598/Documents/combined_df.xlsx"
data = pd.read_excel(file_path)

# 选择特征列和目标变量列
X = data[['pident', 'evalue', 'bitscore', 'distance', 'Alignment_Score']]
y = data['PRF']  # 二分类标签

# 数据标准化
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

# 创建SVM模型，使用RBF核，设置probability=True以输出概率
model = SVC(kernel='rbf', C=1.0, gamma='scale', probability=True)

# 创建分层交叉验证对象
stratified_kfold = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

# 初始化列表存储每个fold的AUC和PR AUC
roc_auc_scores = []
pr_auc_scores = []

for train_index, test_index in stratified_kfold.split(X_scaled, y):
    X_train, X_test = X_scaled[train_index], X_scaled[test_index]
    y_train, y_test = y[train_index], y[test_index]

    # 拟合模型
    model.fit(X_train, y_train)

    # 预测概率
    y_prob = model.predict_proba(X_test)[:, 1]

    # 计算AUC-ROC
    roc_auc = roc_auc_score(y_test, y_prob)
    roc_auc_scores.append(roc_auc)

    # 计算Precision-Recall曲线和AUC
    precision, recall, _ = precision_recall_curve(y_test, y_prob)
    pr_auc = auc(recall, precision)
    pr_auc_scores.append(pr_auc)

    # 绘制ROC曲线
    RocCurveDisplay.from_predictions(y_test, y_prob)
    plt.title(f"ROC Curve for fold with AUC: {roc_auc:.2f}")
    plt.show()

    # 绘制Precision-Recall曲线
    PrecisionRecallDisplay(precision=precision, recall=recall).plot()
    plt.title(f"Precision-Recall Curve for fold with AUC: {pr_auc:.2f}")
    plt.show()

# 打印平均AUC-ROC和PR AUC分数
print(f"Mean AUC-ROC: {np.mean(roc_auc_scores):.2f}")
print(f"Mean Precision-Recall AUC: {np.mean(pr_auc_scores):.2f}")