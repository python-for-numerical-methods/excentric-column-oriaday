import numpy as np
from scipy.optimize import bisect

def find_critical_load(L, E, A, r, c, e, sigma_allow):
    """
    L: אורך במ"מ
    E: מודול אלסטיות ב-MPa
    A: שטח חתך בממ"ר
    r: רדיוס אינרציה במ"מ
    c: מרחק לסיב קיצוני במ"מ
    e: אקסצנטריות במ"מ
    sigma_allow: מאמץ מותר ב-MPa
    
    Return: העומס P בניוטון (float)
    """
    
    # חישוב עומס הקריסה האידיאלי של אוילר (החסם העליון התיאורטי)
    # P_euler = (pi^2 * E * I) / L^2, וכיוון ש- I = A * r^2:
    P_euler = (np.pi**2 * E * A) / (L / r)**2
    
    # פונקציית העזר שאנו רוצים למצוא את השורש שלה: f(P) = 0
    def f(P):
        if P == 0:
            return -sigma_allow  # ב-P=0 המאמץ הוא אפס, ולכן נחזיר את ההפרש
            
        # חישוב הביטוי בתוך הקוסינוס (ברדיאנים)
        theta = (L / (2 * r)) * np.sqrt(P / (E * A))
        
        # חישוב ה-secant (אחד חלקי קוסינוס)
        sec_val = 1.0 / np.cos(theta)
        
        # חישוב המאמץ המקסימלי לפי נוסחת הסקנט
        sigma_max = (P / A) * (1 + (e * c / r**2) * sec_val)
        
        # אנו מחפשים את הנקודה שבה ההפרש שווה לאפס
        return sigma_max - sigma_allow

    # שימוש בשיטת החצייה כדי למצוא את העומס P
    # אנו מחפשים בקטע שבין 0 ל- 99.99% מעומס אוילר (כדי להימנע מחלוקה באפס או אסימפטוטה)
    P_critical = bisect(f, 0, 0.9999 * P_euler)
    
    return float(P_critical)
