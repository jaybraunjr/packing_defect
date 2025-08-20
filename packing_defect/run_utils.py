# packing_defect/run_utils.py

def run_analysis(analyzer):
    """
    Generic runner for any BaseDefectAnalyzer.
    """
    analyzer.run()
    analyzer.plot()
    analyzer.save_results( analyzer.results)
