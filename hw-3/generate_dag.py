import luigi
from graphviz import Digraph
from luigi_pipeline import RunVariantCalling

def generate_dag_image(task_class, filename='dag_output'):
    dot = Digraph(comment='Luigi Pipeline')
    dot.attr('node', shape='box', style='rounded,filled', fillcolor='lightblue')
    
    edges_added = set()

    def add_tasks(task):
        task_id = task.task_family
        dot.node(task_id, task_id)
        
        reqs = luigi.task.flatten(task.requires())
        for req in reqs:
            req_id = req.task_family
            dot.node(req_id, req_id)
            
            edge = (req_id, task_id)
            if edge not in edges_added:
                dot.edge(req_id, task_id)
                edges_added.add(edge)
                
            add_tasks(req)
            
    top_task = task_class()
    add_tasks(top_task)
    
    dot.render(filename, format='png', cleanup=True)
    print(f"DAG generated and saved as {filename}.png")

if __name__ == "__main__":
    generate_dag_image(RunVariantCalling)
