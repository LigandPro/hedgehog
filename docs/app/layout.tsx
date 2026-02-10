import { Footer, Layout, Navbar } from 'nextra-theme-docs'
import { Head } from 'nextra/components'
import { getPageMap } from 'nextra/page-map'
import 'nextra-theme-docs/style.css'

export const metadata = {
  title: {
    default: 'HEDGEHOG',
    template: '%s — HEDGEHOG'
  },
  description: 'Hierarchical Evaluation of Drug GEnerators tHrOugh riGorous filtration'
}

const navbar = (
  <Navbar
    logo={<b>HEDGEHOG</b>}
    projectLink="https://github.com/LigandPro/hedgehog"
  />
)
const footer = <Footer>MIT {new Date().getFullYear()} © Ligand Pro</Footer>

export default async function RootLayout({ children }: { children: React.ReactNode }) {
  return (
    <html lang="en" dir="ltr" suppressHydrationWarning>
      <Head />
      <body>
        <Layout
          navbar={navbar}
          pageMap={await getPageMap()}
          docsRepositoryBase="https://github.com/LigandPro/hedgehog/tree/main/docs"
          footer={footer}
          editLink="Edit this page on GitHub"
        >
          {children}
        </Layout>
      </body>
    </html>
  )
}
